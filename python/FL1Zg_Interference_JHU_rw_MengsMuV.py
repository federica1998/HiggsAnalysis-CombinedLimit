from HiggsAnalysis.CombinedLimit.PhysicsModel import *
import re

class FL1Zg_Interference_JHU_rw_MengsMuV(PhysicsModel):
    def __init__(self):
        self.altSignal = "ALT_L1Zg"

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
            "sigmaa1a3int_VBF": 0,
            "sigmaa1a3int_ZH": 0,
            "sigmaa1a3int_WH": 0,
            "sigmaa1a2int_VBF": 2207.73 ,
            "sigmaa1a2int_ZH": 4258.966,
            "sigmaa1a2int_WH": 16486.68,
            "sigmaa1L1int_VBF": 2861.21349769,
            "sigmaa1L1int_ZH": 6852.307,
            "sigmaa1L1int_WH": 25302.37,
	    "sigmaL1Zg_HZZ": 5.013e-06,
	    "sigmaL1Zg_VBF": 4.9845301e-05,
	    "sigmaL1Zg_ZH": 0.021825343267033705,
	    "sigmaL1Zg_WH": 0,
	    "sigmaa1L1Zgint_VBF": 1410.5494,
	    "sigmaa1L1Zgint_ZH": 21952.490480117736,
	    "sigmaa1L1Zgint_WH": 30998.54,  
        }

        self.modelBuilder.doVar("CMS_zz4l_fai1[0.0,-1.0,1.0]");
#        self.modelBuilder.doVar("fa3_ggH[0.0,0.0,1.0]");
        self.modelBuilder.doVar("muf[1.0,0,10]");
        self.modelBuilder.doVar("muV[1.0,0.0,10.0]");
        self.modelBuilder.doSet("POI","CMS_zz4l_fai1,muV,muf")
#        self.modelBuilder.out.var("muf").setAttribute("flatParam")
        
        self.modelBuilder.doVar('expr::a1("sqrt(1-abs(@0))", CMS_zz4l_fai1)')
        self.modelBuilder.doVar('expr::a3("(@0>0 ? 1 : -1) * sqrt(abs(@0)*{sigma1_HZZ}/{sigmaL1Zg_HZZ})", CMS_zz4l_fai1)'.format(**xsecs))

        self.modelBuilder.factory_('expr::muVc("@1/(1+40.*abs(@0))", CMS_zz4l_fai1,muV)'.format(**xsecs))	
        self.modelBuilder.factory_('expr::smCoupling_VBF("@0*@1**2 - @0*@1*@2*sqrt({sigmaL1Zg_VBF}/{sigma1_VBF})", muVc,a1,a3)'.format(**xsecs))
        self.modelBuilder.factory_('expr::smCoupling_ZH("@0*@1**2 - @0*@1*@2*sqrt({sigmaL1Zg_ZH}/{sigma1_ZH})", muVc,a1,a3)'.format(**xsecs))
        self.modelBuilder.factory_('expr::smCoupling_WH("@0*@1**2 - @0*@1*@2*sqrt({sigmaL1Zg_WH}/{sigma1_WH})", muVc,a1,a3)'.format(**xsecs))


        self.modelBuilder.factory_('expr::bsmCoupling_VBF("@0*@1**2*{sigmaL1Zg_VBF}/{sigma1_VBF} - @0*@1*@2*sqrt({sigmaL1Zg_VBF}/{sigma1_VBF})", muVc,a3,a1)'.format(**xsecs))
        self.modelBuilder.factory_('expr::bsmCoupling_ZH("@0*@1**2*{sigmaL1Zg_ZH}/{sigma1_ZH} - @0*@1*@2*sqrt({sigmaL1Zg_ZH}/{sigma1_ZH})", muVc,a3,a1)'.format(**xsecs))
        self.modelBuilder.factory_('expr::bsmCoupling_WH("@0*@1**2*{sigmaL1Zg_WH}/{sigma1_WH} - @0*@1*@2*sqrt({sigmaL1Zg_WH}/{sigma1_WH})", muVc,a3,a1)'.format(**xsecs))


        self.modelBuilder.factory_('expr::intCoupling_VBF("@0*@1*@2*sqrt({sigmaL1Zg_VBF}/{sigma1_VBF})*{sigmaa1L1Zgint_VBF}/{sigma1_VBF}", muVc,a1,a3)'.format(**xsecs))
        self.modelBuilder.factory_('expr::intCoupling_ZH("@0*@1*@2*sqrt({sigmaL1Zg_ZH}/{sigma1_ZH})*{sigmaa1L1Zgint_ZH}/{sigma1_ZH}", muVc,a1,a3)'.format(**xsecs))
        self.modelBuilder.factory_('expr::intCoupling_WH("@0*@1*@2*sqrt({sigmaL1Zg_WH}/{sigma1_WH})*{sigmaa1L1Zgint_WH}/{sigma1_WH}", muVc,a1,a3)'.format(**xsecs))


    def getYieldScale(self,bin,process):

        if not self.DC.isSignal[process]: return 1

        years = ["2016preVFP","2016postVFP","2017","2018"]
        scale = '1'
        if process in ["qqH_%s_hgg"%y for y in years]:
            scale = 'smCoupling_VBF'
        if process in ["w_minus_%s_hgg"%y for y in years]:
            scale = 'smCoupling_WH'
        if process in ["w_plus_%s_hgg"%y for y in years]:
            scale = 'smCoupling_WH'
        if process in ["zh_%s_hgg"%y for y in years]:
            scale = 'smCoupling_ZH'
        if process in ["qqH_%s_%s_hgg"%(self.altSignal,y) for y in years]:
            scale = 'bsmCoupling_VBF'
        if process in ["wh_%s_%s_hgg"%(self.altSignal,y) for y in years]:
            scale = 'bsmCoupling_WH'
        if process in ["zh_%s_%s_hgg"%(self.altSignal,y) for y in years]:
            scale = 'bsmCoupling_ZH'
        if process in ["qqH_%sf05_%s_hgg"%(self.altSignal,y) for y in years]:
            scale = 'intCoupling_VBF'
        if process in ["zh_%sf05_%s_hgg"%(self.altSignal,y) for y in years]:
            scale = 'intCoupling_ZH'
        if process in ["wh_%sf05_%s_hgg"%(self.altSignal,y) for y in years]:
            scale = 'intCoupling_WH'
        ##############################
        if process in ["ttH_%s_hgg"%y for y in years]:
            scale = 'muf'
        print "Bin/Process %s/%s will get scaled by %s"%(bin,process,scale)
        return scale

FL1Zg_Interference_JHU_rw_MengsMuV = FL1Zg_Interference_JHU_rw_MengsMuV()
