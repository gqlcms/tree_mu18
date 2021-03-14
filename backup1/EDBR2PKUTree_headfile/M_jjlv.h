#ifndef _Mjjlv_
#define _Mjjlv_

void EDBR2PKUTree::Mjjlv() {
        m_JJlv             = Double_t(candMasspuppiJEC);
        m_Jlv_new=-99;
        m_JJlv_new=-99;
        m_Jlv_JEC_up=-99;
        m_Jlv_JEC_down=-99;
        m_Jlv_JER_up=-99;
        m_Jlv_JER_down=-99;
        m_JJlv_JEC_up=-99;
        m_JJlv_JEC_down=-99;
        m_JJlv_JER_up=-99;
        m_JJlv_JER_down=-99;
        
        if (candMasspuppiJEC_new>0) m_JJlv_new             = Double_t(candMasspuppiJEC_new);
        if (candMasspuppiJEC_JEC_up>0) m_JJlv_JEC_up             = Double_t(candMasspuppiJEC_JEC_up);
        if (candMasspuppiJEC_JEC_down>0) m_JJlv_JEC_down             = Double_t(candMasspuppiJEC_JEC_down);
        if (candMasspuppiJEC_JER_up>0) m_JJlv_JER_up             = Double_t(candMasspuppiJEC_JER_up);
        if (candMasspuppiJEC_JER_down>0) m_JJlv_JER_down             = Double_t(candMasspuppiJEC_JER_down);
        
        
        m_Jlv = m_jlv;
        if (m_jlv_new>0) m_Jlv_new = m_jlv_new;
        if (m_jlv_JEC_up>0) m_Jlv_JEC_up = m_jlv_JEC_up;
        if (m_jlv_JEC_down>0) m_Jlv_JEC_down = m_jlv_JEC_down;
        if (m_jlv_JER_up>0) m_Jlv_JER_up = m_jlv_JER_up;
        if (m_jlv_JER_down>0) m_Jlv_JER_down = m_jlv_JER_down;
        
        
}

#endif