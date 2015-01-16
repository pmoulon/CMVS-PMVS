#include <iostream>
#include <fstream>
#include "../image/image.h"
#include "detectFeatures.h"
#include "harris.h"
#include "dog.h"
#include "point.h"

using namespace PMVS3;
using namespace std;
using namespace Image;

CdetectFeatures::CdetectFeatures(void) {
  mtx_init(&m_rwlock, mtx_plain | mtx_recursive);
}

CdetectFeatures::~CdetectFeatures() {
  mtx_destroy(&m_rwlock);
}

void CdetectFeatures::run(const CphotoSetS& pss, const int num,
                          const int csize, const int level,
                          const int CPU) {
  m_ppss = &pss;
  m_csize = csize;
  m_level = level;
  m_CPU = CPU;

  m_points.clear();
  m_points.resize(num);
  
  //----------------------------------------------------------------------
  for (int index = 0; index < num; ++index)
    m_jobs.push_back(index);
  
  vector<thrd_t> threads(m_CPU);
  for (int i = 0; i < m_CPU; ++i)
    thrd_create(&threads[i], &runThreadTmp, (void*)this);
  for (int i = 0; i < m_CPU; ++i)
    thrd_join(threads[i], NULL);
  //----------------------------------------------------------------------
  cerr << "done" << endl;
}

int CdetectFeatures::runThreadTmp(void* arg) {
  CdetectFeatures* detectFeatures = (CdetectFeatures*)arg;  
  detectFeatures->runThread();
  return 0;
}

void CdetectFeatures::runThread(void) {
  while (1) {
    int index = -1;
    mtx_lock(&m_rwlock);
    if (!m_jobs.empty()) {
      index = m_jobs.front();
      m_jobs.pop_front();
    }
    mtx_unlock(&m_rwlock);
    if (index == -1)
      break;
    
    const int image = m_ppss->m_images[index];
    cerr << image << ' ' << flush;

    //?????????????  May need file lock, because targetting images
    //should not overlap among multiple processors.    
    char buffer[1024];
    sprintf(buffer, "%smodels/%08d.affin%d", m_ppss->m_prefix.c_str(), image, m_level);
    ifstream ifstr;
    ifstr.open(buffer);
    if (ifstr.is_open()) {
      ifstr.close();
      continue;
    }
    ifstr.close();
    
    //----------------------------------------------------------------------
    // parameters
    // for harris
    const float sigma = 4.0f;
    // for DoG
    const float firstScale = 1.0f;    const float lastScale = 3.0f;

    //----------------------------------------------------------------------
    // Harris
    {
      Charris harris;
      multiset<Cpoint> result;
      harris.run(m_ppss->m_photos[index].getImage(m_level),
                 m_ppss->m_photos[index].Cimage::getMask(m_level),
                 m_ppss->m_photos[index].Cimage::getEdge(m_level),
                 m_ppss->m_photos[index].getWidth(m_level),
                 m_ppss->m_photos[index].getHeight(m_level), m_csize, sigma, result);
      
      multiset<Cpoint>::reverse_iterator rbegin = result.rbegin();
      while (rbegin != result.rend()) {
        m_points[index].push_back(*rbegin);
        rbegin++;
      }
    }

    //----------------------------------------------------------------------
    // DoG
    {
      Cdog dog;
      multiset<Cpoint> result;
      dog.run(m_ppss->m_photos[index].getImage(m_level),
              m_ppss->m_photos[index].Cimage::getMask(m_level),
              m_ppss->m_photos[index].Cimage::getEdge(m_level),
              m_ppss->m_photos[index].getWidth(m_level),
              m_ppss->m_photos[index].getHeight(m_level),
              m_csize, firstScale, lastScale, result);
      
      multiset<Cpoint>::reverse_iterator rbegin = result.rbegin();      
      while (rbegin != result.rend()) {
        m_points[index].push_back(*rbegin);
        rbegin++;
      }
    }
  }
}
