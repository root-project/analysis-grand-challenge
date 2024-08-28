import ROOT

def MuteTool():
    ROOT.RooMsgService.instance().getStream(1).removeTopic(ROOT.RooFit.Minimization);
    ROOT.RooMsgService.instance().getStream(1).removeTopic(ROOT.RooFit.NumIntegration);
    ROOT.RooMsgService.instance().getStream(1).removeTopic(ROOT.RooFit.Eval);


