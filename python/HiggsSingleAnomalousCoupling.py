from HiggsAnalysis.CombinedLimit.PhysicsModel import *
import re

class FA3_Interference_JHU_ggHSyst_rw_MengsMuV_HeshyXsec_ggHInt_ggHphase(PhysicsModel):
    def __init__(self):
        self.altSignal = "ALT0M"

    def setPhysicsOptions(self,physOptions):
        for po in physOptions:
            if po.startswith("altSignal="): self.altSignal = po.split("=")[1]

    def doParametersOfInterest(self):
        """Create POI and other parameters, and define the POI set."""
        xsecs = {
            "sigma1_HZZ": 290.58626,
            #          "sigma3_HZZ": 581.17253,
            "sigma3_HZZ": 44.670158,
            "sigma1_VBF": 968.674,
            "sigma3_VBF": 10909.54,
            "sigma1_ZH":  9022.36,
            "sigma3_ZH":  434763.7,
            "sigma1_WH":  30998.54,
            "sigma3_WH":  2028656,
            "sigma2_HZZ": 105.85594,
            "sigmaL1_HZZ": 1.9846071e-06,
            "sigma2_VBF": 13102.71,
            "sigmaL1_VBF": 2.08309E-4,
            "sigma2_ZH": 713123,
            "sigmaL1_ZH": 33652.46e-6,
            "sigma2_WH": 3106339,
            "sigmaL1_WH": 11234.91e-5,
            "sigmaa1a3int_VBF": 1937.15,
            "sigmaa1a3int_ZH": 18044.72,
            "sigmaa1a3int_WH": 61997.07,
            "sigma_SM_ggH": 15980,
            "sigma_BSM_ggH": 15981,
        }

        self.modelBuilder.doVar("CMS_zz4l_fai1[0.0,-1.0,1.0]");
        self.modelBuilder.doVar("fa3_ggH[0.0,-1.0,1.0]");
        self.modelBuilder.doVar("muf[1.0,0,10]");
        self.modelBuilder.doVar("mutau[1.0,0,10]");
        self.modelBuilder.doVar("muV[1.0,0.0,10.0]");
        self.modelBuilder.doSet("POI","CMS_zz4l_fai1,fa3_ggH,muV,muf")
#        self.modelBuilder.out.var("muf").setAttribute("flatParam")
        
        self.modelBuilder.doVar('expr::a1("sqrt(1-abs(@0))", CMS_zz4l_fai1)')
        self.modelBuilder.doVar('expr::a3("(@0>0 ? 1 : -1) * sqrt(abs(@0)*{sigma1_HZZ}/{sigma3_HZZ})", CMS_zz4l_fai1)'.format(**xsecs))

        self.modelBuilder.doVar('expr::a1_ggH("sqrt(1-abs(@0))", fa3_ggH)')
        #self.modelBuilder.doVar('expr::a3_ggH("sqrt(abs(@0))", fa3_ggH)'.format(**xsecs)) # no ggH phase!
	# with ggH phase:
        self.modelBuilder.doVar('expr::a3_ggH("(@0>0 ? -1 : 1) * sqrt(abs(@0))", fa3_ggH)'.format(**xsecs))

        self.modelBuilder.factory_('expr::muVc("@1/(1+300.*abs(@0))", CMS_zz4l_fai1,muV)'.format(**xsecs))
#        self.modelBuilder.factory_('expr::muVc("@2/(@0**2+@1**2*{sigma3_VBF}/{sigma1_VBF})", a1,a3,muV)'.format(**xsecs))
#        self.modelBuilder.factory_('expr::muVc_ZZ("@2/(@0**4+@1**4*{sigma3_VBF_ZZ}/{sigma1_VBF_ZZ})", a1,a3,muV)'.format(**xsecs))

        self.modelBuilder.factory_('expr::smCoupling_VBF("@0*@1**2 - @0*@1*@2*sqrt({sigma3_VBF}/{sigma1_VBF})", muVc,a1,a3)'.format(**xsecs))
        self.modelBuilder.factory_('expr::smCoupling_ZH("@0*@1**2 - @0*@1*@2*sqrt({sigma3_ZH}/{sigma1_ZH})", muVc,a1,a3)'.format(**xsecs))
        self.modelBuilder.factory_('expr::smCoupling_WH("@0*@1**2 - @0*@1*@2*sqrt({sigma3_WH}/{sigma1_WH})", muVc,a1,a3)'.format(**xsecs))

        self.modelBuilder.factory_('expr::smCoupling_ggH("@0*@1**2 -  @0*@1*@2*sqrt({sigma_BSM_ggH}/{sigma_SM_ggH})", muf,a1_ggH,a3_ggH)'.format(**xsecs))
        self.modelBuilder.factory_('expr::bsmCoupling_ggH("@0*@1**2*{sigma_BSM_ggH}/{sigma_SM_ggH} - @0*@1*@2*sqrt({sigma_BSM_ggH}/{sigma_SM_ggH})", muf,a3_ggH,a1_ggH)'.format(**xsecs))
        self.modelBuilder.factory_('expr::intCoupling_ggH("@0*@1*@2*sqrt({sigma_BSM_ggH}/{sigma_SM_ggH})*2.", muf,a1_ggH,a3_ggH)'.format(**xsecs))

        self.modelBuilder.factory_('expr::smCoupling_ggH_ZZ("@0*@1**2*@2", muf,a1_ggH,mutau)'.format(**xsecs))
        self.modelBuilder.factory_('expr::bsmCoupling_ggH_ZZ("@0*@1**2*@2", muf,a3_ggH,mutau)'.format(**xsecs))

        self.modelBuilder.factory_('expr::smCoupling_VBF_ZZ("@0*@1**2*@2", muVc,a1,mutau)'.format(**xsecs))
        self.modelBuilder.factory_('expr::bsmCoupling_VBF_ZZ("@0*@1**2*@2", muVc,a3,mutau)'.format(**xsecs))

        self.modelBuilder.factory_('expr::bsmCoupling_VBF("@0*@1**2*{sigma3_VBF}/{sigma1_VBF} - @0*@1*@2*sqrt({sigma3_VBF}/{sigma1_VBF})", muVc,a3,a1)'.format(**xsecs))
        self.modelBuilder.factory_('expr::bsmCoupling_ZH("@0*@1**2*{sigma3_ZH}/{sigma1_ZH} - @0*@1*@2*sqrt({sigma3_ZH}/{sigma1_ZH})", muVc,a3,a1)'.format(**xsecs))
        self.modelBuilder.factory_('expr::bsmCoupling_WH("@0*@1**2*{sigma3_WH}/{sigma1_WH} - @0*@1*@2*sqrt({sigma3_WH}/{sigma1_WH})", muVc,a3,a1)'.format(**xsecs))



        self.modelBuilder.factory_('expr::intCoupling_VBF("@0*@1*@2*sqrt({sigma3_VBF}/{sigma1_VBF})*{sigmaa1a3int_VBF}/{sigma1_VBF}", muVc,a1,a3)'.format(**xsecs))
        self.modelBuilder.factory_('expr::intCoupling_ZH("@0*@1*@2*sqrt({sigma3_ZH}/{sigma1_ZH})*{sigmaa1a3int_ZH}/{sigma1_ZH}", muVc,a1,a3)'.format(**xsecs))
        self.modelBuilder.factory_('expr::intCoupling_WH("@0*@1*@2*sqrt({sigma3_WH}/{sigma1_WH})*{sigmaa1a3int_WH}/{sigma1_WH}", muVc,a1,a3)'.format(**xsecs))


    def getYieldScale(self,bin,process):

        if not self.DC.isSignal[process]: return 1

        years = [2016,2017,2018]
        scale = '1'
        if process in ["ggH_%d_hgg"%y for y in years]:
            scale = 'smCoupling_ggH'
        if process in ["ggH_%s_%d_hgg"%(self.altSignal,y) for y in years]: # not yet considered
            scale = 'bsmCoupling_ggH'
        if process in ["qqH_%d_hgg"%y for y in years]:
            scale = 'smCoupling_VBF'
        if process in ["vH_%d_hgg"%y for y in years]: # should separate W and Z
            scale = 'smCoupling_WH'
        if process in ["zH_%d_hgg"%y for y in years]: # not yet there
            scale = 'smCoupling_ZH'
        if process in ["qqH_%s_%d_hgg"%(self.altSignal,y) for y in years]:
            scale = 'bsmCoupling_VBF'
        if process in ["wH_%s_%d_hgg"%(self.altSignal,y) for y in years]: # not yet considered
            scale = 'bsmCoupling_WH'
        if process in ["wH_%s_%d_hgg"%(self.altSignal,y) for y in years]: # not yet considered
            scale = 'bsmCoupling_ZH'
        if process in ["qqH_%sf05_%d_hgg"%(self.altSignal,y) for y in years]:
            scale = 'intCoupling_VBF'
        if process in ["zH_%sf05_%d_hgg"%(self.altSignal,y) for y in years]: # not yet considered
            scale = 'intCoupling_ZH'
        if process in ["wH_%sf05_%d_hgg"%(self.altSignal,y) for y in years]: # not yet considered
            scale = 'intCoupling_WH'
        if process in ["ggH_%sf05_%d_hgg"%(self.altSignal,y) for y in years]: # not yet considered
            scale = 'intCoupling_ggH'
        print "Bin/Process %s/%s will get scaled by %s"%(bin,process,scale)
        return scale


FA3_Interference_JHU_ggHSyst_rw_MengsMuV_HeshyXsec_ggHInt_ggHphase = FA3_Interference_JHU_ggHSyst_rw_MengsMuV_HeshyXsec_ggHInt_ggHphase()
