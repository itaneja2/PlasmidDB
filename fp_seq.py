

fp_seq_map = {} 

#sequence derived from https://www.fpbase.org/protein/mmaple3/
mmaple3_seq = 'MVSKGEETIM SVIKPDMKIK LRMEGNVNGH AFVIEGEGSG KPFEGIQTID LEVKEGAPLP FAYDILTTAF HYGNRVFTKY PRKIPDYFKQ SFPEGYSWER SMTYEDGGIC NATNDITMEE DSFINKIHFK GTNFPPNGPV MQKRTVGWEV STEKMYVRDG VLKGDVKMKL LLKGGSHYRC DFRTTYKVKQ KAVKLPKAHF VDHRIEILSH DKDYNKVKLY EHAVARNSTD SMDELYK'.replace(' ','')
fp_seq_map['mMaple3'] = mmaple3_seq

mcherry_seq = 'MVSKGEEDNM AIIKEFMRFK VHMEGSVNGH EFEIEGEGEG RPYEGTQTAK LKVTKGGPLP FAWDILSPQF MYGSKAYVKH PADIPDYLKL SFPEGFKWER VMNFEDGGVV TVTQDSSLQD GEFIYKVKLR GTNFPSDGPV MQKKTMGWEA SSERMYPEDG ALKGEIKQRL KLKDGGHYDA EVKTTYKAKK PVQLPGAYNV NIKLDITSHN EDYTIVEQYE RAEGRHSTGG MDELYK'.replace(' ','')
fp_seq_map['mCherry'] = mcherry_seq

eyfp_seq = 'MVSKGEELFT GVVPILVELD GDVNGHKFSV SGEGEGDATY GKLTLKFICT TGKLPVPWPT LVTTFGYGLQ CFARYPDHMK QHDFFKSAMP EGYVQERTIF FKDDGNYKTR AEVKFEGDTL VNRIELKGID FKEDGNILGH KLEYNYNSHN VYIMADKQKN GIKVNFKIRH NIEDGSVQLA DHYQQNTPIG DGPVLLPDNH YLSYQSALSK DPNEKRDHMV LLEFVTAAGI TLGMDELYK'.replace(' ','')
fp_seq_map['eYFP'] = eyfp_seq

sfgfp_seq = 'MSKGEELFTG VVPILVELDG DVNGHKFSVR GEGEGDATNG KLTLKFICTT GKLPVPWPTL VTTLTYGVQC FSRYPDHMKR HDFFKSAMPE GYVQERTISF KDDGTYKTRA EVKFEGDTLV NRIELKGIDF KEDGNILGHK LEYNFNSHNV YITADKQKNG IKANFKIRHN VEDGSVQLAD HYQQNTPIGD GPVLLPDNHY LSTQSVLSKD PNEKRDHMVL LEFVTAAGIT HGMDELYK'.replace(' ','')
fp_seq_map['sfGFP'] = sfgfp_seq

redfast_seq = 'MEHVAFGSEDIENTLAKMDDGQLDGLALGAIQLDGDGNILQYNAAQGDITGADPKQVIGKNFFKDVAPGTDSPEFYGKFKVGVASGNLNTMFEWMIPTNRGPTKVKVHMKKALSGDSYWVFVKRV'.replace(' ','')
fp_seq_map['redFAST'] = redfast_seq 

mscarlet_seq = 'MVSKGEAVIK EFMRFKVHME GSMNGHEFEI EGEGEGRPYE GTQTAKLKVT KGGPLPFSWD ILSPQFMYGS RAFTKHPADI PDYYKQSFPE GFKWERVMNF EDGGAVTVTQ DTSLEDGTLI YKVKLRGTNF PPDGPVMQKK TMGWEASTER LYPEDGVLKG DIKMALRLKD GGRYLADFKT TYKAKKPVQM PGAYNVDRKL DITSHNEDYT VVEQYERSEG RHSTGGMDEL YK'.replace(' ','')
fp_seq_map['mScarlet'] = mscarlet_seq

meos_seq = 'MSAIKPDMKI NLRMEGNVNG HHFVIDGDGT GKPFEGKQSM DLEVKEGGPL PFAFDILTTA FHYGNRVFAE YPDHIQDYFK QSFPKGYSWE RSLTFEDGGI CIARNDITME GDTFYNKVRF HGTNFPANGP VMQKKTLKWE PSTEKMYVRD GVLTGDIHMA LLLEGNAHYR CDFRTTYKAK EKGVKLPGYH FVDHCIEILS HDKDYNKVKL YEHAVAHSGL PDNARR'.replace(' ','')
fp_seq_map['mEos'] = meos_seq

meos32_seq = 'MSAIKPDMKI KLRMEGNVNG HHFVIDGDGT GKPFEGKQSM DLEVKEGGPL PFAFDILTTA FHYGNRVFAK YPDNIQDYFK QSFPKGYSWE RSLTFEDGGI CNARNDITME GDTFYNKVRF YGTNFPANGP VMQKKTLKWE PSTEKMYVRD GVLTGDIEMA LLLEGNAHYR CDFRTTYKAK EKGVKLPGAH FVDHCIEILS HDKDYNKVKL YEHAVAHSGL PDNARR'.replace(' ','')
fp_seq_map['mEos3.2'] = meos32_seq

mscarlet_seq = 'MVSKGEAVIK EFMRFKVHME GSMNGHEFEI EGEGEGRPYE GTQTAKLKVT KGGPLPFSWD ILSPQFMYGS RAFTKHPADI PDYYKQSFPE GFKWERVMNF EDGGAVTVTQ DTSLEDGTLI YKVKLRGTNF PPDGPVMQKK TMGWEASTER LYPEDGVLKG DIKMALRLKD GGRYLADFKT TYKAKKPVQM PGAYNVDRKL DITSHNEDYT VVEQYERSEG RHSTGGMDEL YK'.replace(' ','')
fp_seq_map['mScarlet'] = mscarlet_seq

mscarlet3_seq = 'MDSTEAVIKE FMRFKVHMEG SMNGHEFEIE GEGEGRPYEG TQTAKLRVTK GGPLPFSWDI LSPQFMYGSR AFTKHPADIP DYWKQSFPEG FKWERVMNFE DGGAVSVAQD TSLEDGTLIY KVKLRGTNFP PDGPVMQKKT MGWEASTERL YPEDVVLKGD IKMALRLKDG GRYLADFKTT YRAKKPVQMP GAFNIDRKLD ITSHNEDYTV VEQYERSVAR HSTGGSGGS'.replace(' ','')
fp_seq_map['mScarlet3'] = mscarlet3_seq 

mscarleti3_seq = 'MDSTEAVIKE FMRFKVHMEG SMNGHEFEIE GEGEGRPYEG TQTAKLKVTK GGPLPFSWDI LSPQFMYGSR AFIKHPADIP DYWKQSFPEG FKWERVMIFE DGGTVSVTQD TSLEDGTLIY KVKLRGGNFP PDGPVMQKRT MGWEASTERL YPEDVVLKGD IKMALRLKDG GRYLADFKTT YKAKKPVQMP GAFNIDRKLD ITSHNEDYTV VEQYERSVAR HSTGGSGGS'.replace(' ','')
fp_seq_map['mScarlet-I3'] = mscarleti3_seq 

redfast_seq = 'MEHVAFGSEDIENTLAKMDDGQLDGLALGAIQLDGDGNILQYNAAQGDITGADPKQVIGKNFFKDVAPGTDSPEFYGKFKVGVASGNLNTMFEWMIPTNRGPTKVKVHMKKALSGDSYWVFVKRV'
fp_seq_map['redFast'] = redfast_seq

halotag_seq = 'MAEIGTGFPFDPHYVEVLGERMHYVDVGPRDGTPVLFLHGNPTSSYVWRNIIPHVAPTHRCIAPDLIGMGKSDKPDLGYFFDDHVRFMDAFIEALGLEEVVLVIHDWGSALGFHWAKRNPERVKGIAFMEFIRPIPTWDEWPEFARETFQAFRTTDVGRKLIIDQNVFIEGTLPMGVVRPLTEVEMDHYREPFLNPVDREPLWRFPNELPIAGEPANIVALVEEYMDWLHQSPVPKLLFWGTPGVLIPPAEAARLAKSLPNCKAVDIGPGLNLLQEDNPDLIGSEIARWLSTLEISG'
fp_seq_map['HaloTag'] = halotag_seq 

egfp_seq = 'MVSKGEELFT GVVPILVELD GDVNGHKFSV SGEGEGDATY GKLTLKFICT TGKLPVPWPT LVTTLTYGVQ CFSRYPDHMK QHDFFKSAMP EGYVQERTIF FKDDGNYKTR AEVKFEGDTL VNRIELKGID FKEDGNILGH KLEYNYNSHN VYIMADKQKN GIKVNFKIRH NIEDGSVQLA DHYQQNTPIG DGPVLLPDNH YLSTQSALSK DPNEKRDHMV LLEFVTAAGI TLGMDELYK'.replace(' ','')
fp_seq_map['EGFP'] = egfp_seq 

megfp_seq = 'MVSKGEELFT GVVPILVELD GDVNGHKFSV SGEGEGDATY GKLTLKFICT TGKLPVPWPT LVTTLTYGVQ CFSRYPDHMK QHDFFKSAMP EGYVQERTIF FKDDGNYKTR AEVKFEGDTL VNRIELKGID FKEDGNILGH KLEYNYNSHN VYIMADKQKN GIKVNFKIRH NIEDGSVQLA DHYQQNTPIG DGPVLLPDNH YLSTQSKLSK DPNEKRDHMV LLEFVTAAGI TLGMDELYK'.replace(' ','')
fp_seq_map['mEGFP'] = megfp_seq 

mkeima_seq = 'MVSVIAKQMT YKVYMSGTVN GHYFEVEGDG KGKPYEGEQT VKLTVTKGGP LPFAWDILSP QLQYGSIPFT KYPEDIPDYF KQSFPEGYTW ERSMNFEDGA VCTVSNDSSI QGNCFIYNVK ISGENFPPNG PVMQKKTQGW EPSTERLFAR DGMLIGNDYM ALKLEGGGHY LCEFKSTYKA KKPVRMPGRH EIDRKLDVTS HNRDYTSVEQ CEIAIARHSL LG'.replace(' ','')
fp_seq_map['mKeima'] = mkeima_seq

sapphire_seq = 'MSKGEELFTG VVPILVELDG DVNGHKFSVS GEGEGDATYG KLTLKFICTT GKLPVPWPTL VTTFSYGVQC FARYPDHMKQ HDFFKSAMPE GYVQERTIFF KDDGNYKTRA EVKFEGDTLV NRIELKGIDF KEDGNILGHK LEYNFNSHNV YIMADKQKNG IKVNFKIRHN IEDGSVQLAD HYQQNTPIGD GPVLLPDNHY LSIQSALSKD PNEKRDHMVL LEFVTAAGIT LGMDELYK'.replace(' ','')
fp_seq_map['Sapphire'] = sapphire_seq

tsapphire_seq = 'MVSKGEELFT GVVPILVELD GDVNGHKFSV SGEGEGDATY GKLTLKFICT TGKLPVPWPT LVTTFSYGVM VFARYPDHMK QHDFFKSAMP EGYVQERTIF FKDDGNYKTR AEVKFEGDTL VNRIELKGID FKEDGNILGH KLEYNFNSHN VYIMADKQKN GIKANFKIRH NIEDGGVQLA DHYQQNTPIG DGPVLLPDNH YLSIQSALSK DPNEKRDHMV LLEFVTAAGI TLGMDELYK'.replace(' ','')
fp_seq_map['T-Sapphire'] = tsapphire_seq

miRFP703_seq = 'MVAGHASGSP AFGTASHSNC EHEEIHLAGS IQPHGALLVV SEHDHRVIQA SANAAEFLNL GSVLGVPLAE IDGDLLIKIL PHLDPTAEGM PVAVRCRIGN PSTEYCGLMH RPPEGGLIIE LERAGPSIDL SGTLAPALER IRTAGSLRAL CDDTVLLFQQ CTGYDRVMVY RFDEQGHGLV FSECHVPGLE SYFGNRYPSS LVPQMARQLY VRQRVRVLVD VTYQPVPLEP RLSPLTGRDL DMSGCFLRSM SPIHLQFLKD MGVRATLAVS LVVGGKLWGL VVCHHYLPRF IRFELRAICK RLAERIATRI TALES'.replace(' ','')
fp_seq_map['miRFP703'] = miRFP703_seq

mEosEM_seq = 'MVSAIKPDMR IKLRMEGNVN GHHFVIDGEG TGKPYEGKQT MDLEVKEGGP LPFAFDILTT AFHYGNRVFV KYPDNIQDYF KQSFPKGYSW ERSMTFEDGG ICNARNDITM EGDTFYNKVR FYGTNFPANG PVMQKKTLKW EPSTEKMYVR DGVLTGDIEM ALLLEGGAHY RCDFRTTYKA KEKGVKLPGA HFVDHAIEIL SHDKDYNKVK LYEHAVAHSG LPDNARR'.replace(' ','')
fp_seq_map['mEosEM'] = mEosEM_seq

yEGFP_sep = 'MVSKGEELFT GVVPILVELD GDVNGHKFSV SGEGEGDATY GKLTLKFICT TGKLPVPWPT LVTTFGYGVQ CFARYPDHMK QHDFFKSAMP EGYVQERTIF FKDDGNYKTR AEVKFEGDTL VNRIELKGID FKEDGNILGH KLEYNYNSHN VYIMADKQKN GIKVNFKIRH NIEDGSVQLA DHYQQNTPIG DGPVLLPDNH YLSTQSALSK DPNEKRDHMV LLEFVTAAGI THGMDELYK'.replace(' ','')
fp_seq_map['yEGFP'] = yEGFP_sep

msgfp2_seq = 'MDSTESLFTGVVPILVELDGDVNGHKFSVRGEGEGDATNGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTITFKDDGTYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYITADKQKNGIKANFKIRHNVEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSKLSKDPNEKRDHMVLLEFVTAAGIT'
fp_seq_map['msGFP2'] = msgfp2_seq

snaptag_seq = 'MDKDCEMKRTTLDSPLGKLELSGCEQGLHEIKLLGKGTSAADAVEVPAPAAVLGGPEPLMQATAWLNAYFHQPEAIEEFPVPALHHPVFQQESFTRQVLWKLLKVVKFGEVISYQQLAALAGNPAATAAVKTALSGNPVPILIPCHRVVSSSGAVGGYEGGLAVKEWLLAHEGHRLGKPGLG'
fp_seq_map['SnapTag'] = snaptag_seq
