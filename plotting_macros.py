import ROOT
import os
from ROOT import TColor
from analyzeZX import setCavasAndStyles

def CRPlot():
    c1 = ROOT.TCanvas()
    setCavasAndStyles("c1",c1,"")   
    
    histName = ["2P2F","2P2F_4e","2P2F_4mu","2P2F_2e2mu","2P2F_2mu2e","3P1F","3P1F_4e","3P1F_4mu","3P1F_2e2mu","3P1F_2mu2e"]
 #   Fillcolors = ["#cc0099","#99ccff","#996666","#669966"]
 #   Linecolors = ["#990066","#000099","#5f3f3f","#003300"]
    StackTitle = "Proba"
    #ProcessNames = ["Data","WZ","ZZ","TT","DY50"]
    #LabelNames = ["Data","#WZ#","#Z\\gamma^*,ZZ#","#t\\bar{t}+jets#","#Z+jets#"]

    Fillcolors = ["#99ccff","#cc0099","#996666","#669966"]
    Linecolors = ["#000099","#990066","#5f3f3f","#003300"]

    ProcessNames = ["Data","ZZ","WZ","TT","DY50"]
    LabelNames = ["Data","#Z\\gamma^*,ZZ#","#WZ#","#t\\bar{t}+jets#","#Z+jets#"]

    
    fileSize = 5
    plotSize = 10
    variableN = "ptl3"
    for p in range(plotSize):
        print( "Ulazi u plotting" ) 
        leg = ROOT.TLegend(0.45,0.5,1.05,0.85)
        leg.SetBorderSize(0)
        leg.SetLineColor(1)
        leg.SetLineStyle(1)
        leg.SetLineWidth(1)
        leg.SetFillColor(0)
        leg.SetFillStyle(0)
    
        if (p>4):
            entry =leg.AddEntry("NULL","Control Region "+ histName[p] ,"h")
        else:
            entry = leg.AddEntry("NULL","Control Region "+ histName[p] ,"h")

        entry.SetLineColor(1)
        entry.SetLineStyle(1)
        entry.SetLineWidth(1)
        entry.SetMarkerColor(1)
        entry.SetMarkerStyle(21)
        entry.SetMarkerSize(1)
        entry.SetTextFont(62)
        entry=leg.AddEntry("NULL","Data","LP")
        entry.SetLineColor(1)
        entry.SetLineStyle(1)
        entry.SetLineWidth(1)
        entry.SetMarkerColor(1)
        entry.SetMarkerStyle(20)
        entry.SetMarkerSize(1)
        entry.SetTextFont(42)
        entry=leg.AddEntry("NULL","Z + jets","F")
        
        ci = TColor.GetColor("#669966")
        entry.SetFillColor(ci)
        entry.SetFillStyle(1001)
        
        ci = TColor.GetColor("#003300")
        entry.SetLineColor(ci)
        entry.SetLineStyle(1)
        entry.SetLineWidth(1)
        entry.SetMarkerColor(1)
        entry.SetMarkerStyle(21)
        entry.SetMarkerSize(1)
        entry.SetTextFont(42)
        entry=leg.AddEntry("NULL","t#bar{t} + jets","F")
        
        ci = TColor.GetColor("#996666")
        entry.SetFillColor(ci)
        entry.SetFillStyle(1001)
        
        ci = TColor.GetColor("#5f3f3f")
        entry.SetLineColor(ci)
        entry.SetLineStyle(1)
        entry.SetLineWidth(1)
        entry.SetMarkerColor(1)
        entry.SetMarkerStyle(21)
        entry.SetMarkerSize(1)
        entry.SetTextFont(42)
        entry=leg.AddEntry("NULL","WZ","F")
        
        ci = TColor.GetColor("#cc0099")
        entry.SetFillColor(ci)
        entry.SetFillStyle(1001)
        
        ci = TColor.GetColor("#990066")
        entry.SetLineColor(ci)
        entry.SetLineStyle(1)
        entry.SetLineWidth(1)
        entry.SetMarkerColor(1)
        entry.SetMarkerStyle(21)
        entry.SetMarkerSize(1)
        entry.SetTextFont(42)
        entry=leg.AddEntry("NULL","Z#gamma*,ZZ","F")
        
        ci = TColor.GetColor("#99ccff")
        entry.SetFillColor(ci)
        entry.SetFillStyle(1001)
        
        if (p>4):

            ci = TColor.GetColor("#000099")
            entry.SetLineColor(ci)
            entry.SetLineStyle(1)
            entry.SetLineWidth(1)
            entry.SetMarkerColor(1)
            entry.SetMarkerStyle(21)
            entry.SetMarkerSize(1)
            entry.SetTextFont(42)
            entry=leg.AddEntry("NULL","2P2F contribution","F")
            entry.SetFillStyle(4000)
            entry.SetLineColor(6)
            entry.SetLineStyle(1)
            entry.SetLineWidth(2)
            entry.SetMarkerColor(1)
            entry.SetMarkerStyle(21)
            entry.SetMarkerSize(1)
            entry.SetTextFont(42)        
        
        
        
        A = ROOT.THStack(StackTitle,StackTitle)
        
        Proba_stack_5_stack_1 = ROOT.TH1F("Proba_stack_5_stack_1","Proba",100,0,2000)
        Proba_stack_5_stack_1.SetMinimum(-5.608576)
        Proba_stack_5_stack_1.SetMaximum(51.10072)
        Proba_stack_5_stack_1.SetDirectory(0)
        Proba_stack_5_stack_1.SetStats(0)
        Proba_stack_5_stack_1.SetLineColor(TColor.GetColor("#000099"))
        Proba_stack_5_stack_1.SetLineStyle(0)
        Proba_stack_5_stack_1.SetMarkerStyle(20)
        Proba_stack_5_stack_1.GetXaxis().SetLabelFont(42)
        Proba_stack_5_stack_1.GetXaxis().SetLabelOffset(0.007)
        Proba_stack_5_stack_1.GetXaxis().SetLabelSize(0.05)
        Proba_stack_5_stack_1.GetXaxis().SetTitleSize(0.06)
        Proba_stack_5_stack_1.GetXaxis().SetTitleOffset(0.9)
        Proba_stack_5_stack_1.GetXaxis().SetTitleFont(42)
        Proba_stack_5_stack_1.GetYaxis().SetLabelFont(42)
        Proba_stack_5_stack_1.GetYaxis().SetLabelOffset(0.007)
        Proba_stack_5_stack_1.GetYaxis().SetLabelSize(0.05)
        Proba_stack_5_stack_1.GetYaxis().SetTitleSize(0.06)
        Proba_stack_5_stack_1.GetYaxis().SetTitleOffset(1.25)
        Proba_stack_5_stack_1.GetYaxis().SetTitleFont(42)
        Proba_stack_5_stack_1.GetZaxis().SetLabelFont(42)
        Proba_stack_5_stack_1.GetZaxis().SetLabelOffset(0.007)
        Proba_stack_5_stack_1.GetZaxis().SetLabelSize(0.05)
        Proba_stack_5_stack_1.GetZaxis().SetTitleSize(0.06)
        Proba_stack_5_stack_1.GetZaxis().SetTitleFont(42)
        A.SetHistogram(Proba_stack_5_stack_1)
        
        _file = []
        histPlot = []
        for i in range(fileSize):

            if (i == 0):
                Title = "Data"
            else:
                Title = "MC"

            _file.append(ROOT.TFile("Hist_"+Title+"_"+variableN+"_"+ProcessNames[i]+".root"))
            
            if (i==0) :
                dataPlot = _file[i].Get("h1D_m4l_"+histName[p])
                histPlot.append(dataPlot)

                print( "Data " + histName[p]+ " := "  + str(dataPlot.Integral()) )
                
            
            else:
                histPlot.append(_file[i].Get("h1D_m4l_"+histName[p]))
                histPlot[i].SetFillColor(TColor.GetColor(Fillcolors[i-1]))
                histPlot[i].SetLineColor(TColor.GetColor(Linecolors[i-1]))
                histPlot[i].SetFillStyle(1001)
                
                histPlot[i].SetLineStyle(0)
                histPlot[i].SetMarkerStyle(20)
                histPlot[i].GetXaxis().SetLabelFont(42)
                histPlot[i].GetXaxis().SetLabelOffset(0.007)
                histPlot[i].GetXaxis().SetLabelSize(0.05)
                histPlot[i].GetXaxis().SetTitleSize(0.06)
                histPlot[i].GetXaxis().SetTitleOffset(0.9)
                histPlot[i].GetXaxis().SetTitleFont(42)
                histPlot[i].GetYaxis().SetLabelFont(42)
                histPlot[i].GetYaxis().SetLabelOffset(0.007)
                histPlot[i].GetYaxis().SetLabelSize(0.05)
                histPlot[i].GetYaxis().SetTitleSize(0.06)
                histPlot[i].GetYaxis().SetTitleOffset(1.25)
                histPlot[i].GetYaxis().SetTitleFont(42)
                histPlot[i].GetZaxis().SetLabelFont(42)
                histPlot[i].GetZaxis().SetLabelOffset(0.007)
                histPlot[i].GetZaxis().SetLabelSize(0.05)
                histPlot[i].GetZaxis().SetTitleSize(0.06)
                histPlot[i].GetZaxis().SetTitleFont(42)
                
                print( "MC: "+ histName[p] + " := "+ str( histPlot[i].Integral() ) )
                A.Add(histPlot[i])

        print ("Data status " + str(dataPlot))    
        c1 = ROOT.TCanvas("c1", "myPlots",0,67,600,600)
        ROOT.gStyle.SetOptFit(1)
        ROOT.gStyle.SetOptStat(0)
        ROOT.gStyle.SetOptTitle(0)
        c1.Range(-102.5,-10.38415,847.5,69.4939)
        c1.SetFillColor(0)
        c1.SetBorderMode(0)
        c1.SetBorderSize(2)
        c1.SetTickx(1)
        c1.SetTicky(1)
        c1.SetLeftMargin(0.15)
        c1.SetRightMargin(0.05)
        c1.SetTopMargin(0.05)
        c1.SetBottomMargin(0.13)
        c1.SetFrameFillStyle(0)
        c1.SetFrameBorderMode(0)
        c1.SetFrameFillStyle(0)
        c1.SetFrameBorderMode(0)
        
        dataPlot.SetFillColor(TColor.GetColor("#000000"))
        dataPlot.SetLineColor(TColor.GetColor("#000000"))
        dataPlot.SetFillStyle(0)
        dataPlot.SetMarkerStyle(20)
        dataPlot.GetXaxis().SetTitle("m_{4#font[12]{l}} (GeV)")
        dataPlot.GetXaxis().SetRange(3,40)
        dataPlot.GetXaxis().SetLabelFont(42)
        dataPlot.GetXaxis().SetLabelSize(0.05)
        dataPlot.GetXaxis().SetTitleSize(0.06)
        dataPlot.GetXaxis().SetTitleOffset(0.9)
        dataPlot.GetXaxis().SetTitleFont(42)
        dataPlot.GetYaxis().SetTitle("Events / 20 GeV")
        dataPlot.GetYaxis().SetLabelFont(42)
        dataPlot.GetYaxis().SetLabelSize(0.05)
        dataPlot.GetYaxis().SetTitleSize(0.06)
        dataPlot.GetYaxis().SetTitleOffset(1.25)
        dataPlot.GetYaxis().SetTitleFont(42)
        dataPlot.GetZaxis().SetLabelFont(42)
        dataPlot.GetZaxis().SetLabelSize(0.035)
        dataPlot.GetZaxis().SetTitleSize(0.035)
        dataPlot.GetZaxis().SetTitleFont(42)
        
        

        dataPlot.GetXaxis().SetRangeUser(50,800)
        dataPlot.GetYaxis().SetRangeUser(0, 1.3 * dataPlot.GetMaximum())
        dataPlot.Draw("e1 goff")
        
        A.Draw("hist same")
        dataPlot.Draw("same e1 goff")
        if ( p>4 ):
            _file2p2f = ROOT.TFile("estimateZXData.root", "READ")
            histPlot2p2f = _file2p2f.Get("h1D_m4l_Add_"+histName[p-5])
            
            print( "2P2F contr:= " + str(histPlot2p2f.Integral()) ) 
            print ("Adding "+ ProcessNames[1] + " and "+ ProcessNames[2] +" to the 2P2F Contr in 3P1F")
            histPlot2p2f.Add(histPlot[1])
            histPlot2p2f.Add(histPlot[2])
            histPlot2p2f.SetFillColor(TColor.GetColor("#ffffff"))
            histPlot2p2f.SetLineColor(TColor.GetColor("#ff00ff"))
            histPlot2p2f.SetLineWidth(2)
            histPlot2p2f.SetFillStyle(4000)
            histPlot2p2f.Smooth()
            histPlot2p2f.Draw("hist E1 same goff")

        leg.Draw("goff")
        
        tex = ROOT.TLatex(0.95,0.96,"58.8 fb^{-1} (13 TeV)")
        tex.SetNDC()
        tex.SetTextAlign(31)
        tex.SetTextFont(42)
        tex.SetTextSize(0.03)
        tex.SetLineWidth(2)
        tex.Draw("goff")
        tex = ROOT.TLatex(0.15,0.96,"CMS")
        tex.SetNDC()
        tex.SetTextFont(61)
        tex.SetTextSize(0.0375)
        tex.SetLineWidth(2)
        tex.Draw("goff")
        tex = ROOT.TLatex(0.23,0.96,"Preliminary")
        tex.SetNDC()
        tex.SetTextFont(52)
        tex.SetTextSize(0.0285)
        tex.SetLineWidth(2)
        tex.Draw("goff")
        
        c1.SaveAs("CR_plots/"+histName[p]+"_auto.pdf")
        c1.SaveAs("CR_plots/m"+histName[p]+"_auto.C")
        c1.SaveAs("CR_plots/"+histName[p]+"_auto.root")



def CRMatchPlot():
    c1 = ROOT.TCanvas()
    setCavasAndStyles("c1",c1,"") 
    
    histName = ["2P2F","2P2F4e","2P2F4mu","2P2F2e2mu","2P2F2mu2e","3P1F","3P1F4e","3P1F4mu","3P1F2e2mu", "3P1F2mu2e"]
    dataName = ["2P2F","2P2F_4e","2P2F_4mu","2P2F_2e2mu","2P2F_2mu2e","3P1F","3P1F_4e","3P1F_4mu","3P1F_2e2mu", "3P1F_2mu2e"]
    Fillcolors = ["#AF0096","#91FFFF","#0000ff","#669965","#3CB371"]
    Linecolors = ["#AF0096","#91FFFF","#0000ff","#669965","#3CB371"]
    StackTitle = "Proba"
    ProcessNames = ["Data","WZ","ZZ","TT","DY50"]
    
    fileSize = 5
    plotSize = 10
    variableN = "ptl3"
    Title = "Data"
    _file1 = ROOT.TFile("Hist_"+Title+"_"+variableN+"_"+Title+".root")
    _file = ROOT.TFile("Hist_ID_"+variableN+".root")
    for p in range(plotSize):        
        leg = 0
        leg_xl = 0.45
        leg_xr = 0.95
        leg_yb = 0.55
        leg_yt = 0.85 

        leg = ROOT.TLegend(leg_xl,leg_yb,leg_xr,leg_yt) 
     #   leg = ROOT.TLegend(0.45,0.5,1.05,0.85)
        leg.SetBorderSize(0)
        leg.SetLineColor(1)
        leg.SetLineStyle(1)
        leg.SetLineWidth(1)
        leg.SetFillColor(0)
        leg.SetFillStyle(0)
        A = ROOT.THStack(StackTitle,StackTitle)
        
        dataPlot = _file1.Get("h1D_m4l_"+dataName[p])
        leg.AddEntry(dataPlot,Title,"LP")
        
        histPlotC = _file.Get("h1D_m4l_Conv_"+histName[p])
        histPlotC.SetFillColor(30)
        histPlotC.SetLineColor(1)
        histPlotC.SetFillStyle(1001)
        
        histPlotC.SetLineStyle(0)
        histPlotC.SetMarkerStyle(20)
        
        leg.AddEntry(histPlotC,"Conversions","F")
        A.Add(histPlotC)
        
        histPlotP = _file.Get("h1D_m4l_Prom_"+histName[p])
        histPlotP.SetFillColor(38)
        histPlotP.SetLineColor(1)
        histPlotP.SetFillStyle(1001)
        
        histPlotP.SetLineStyle(0)
        histPlotP.SetMarkerStyle(20)
        leg.AddEntry(histPlotP,"Prompt","F")
        A.Add(histPlotP)
        
        histPlotF = _file.Get("h1D_m4l_Fake_"+histName[p])
        histPlotF.SetFillColor(46)
        histPlotF.SetLineColor(1)
        histPlotF.SetFillStyle(1001)
        
        histPlotF.SetLineStyle(0)
        histPlotF.SetMarkerStyle(20)
        leg.AddEntry(histPlotF,"LF Fakes","F")
        A.Add(histPlotF)
        
        histPlotBD = _file.Get("h1D_m4l_BDfake_"+histName[p])
        histPlotBD.SetFillColor(49)
        histPlotBD.SetLineColor(1)
        histPlotBD.SetFillStyle(1001)
        histPlotBD.SetLineStyle(0)
        histPlotBD.SetMarkerStyle(20)
        leg.AddEntry(histPlotBD,"HF fakes","F")
        A.Add(histPlotBD)
        
        # setup environment & canvas
        
        c1 = ROOT.TCanvas("c1", "myPlots",0,67,600,600)
        ROOT.gStyle.SetOptFit(1)
        ROOT.gStyle.SetOptStat(0)
        ROOT.gStyle.SetOptTitle(0)
        c1.Range(-102.5,-10.38415,847.5,69.4939)
        c1.SetFillColor(0)
        c1.SetBorderMode(0)
        c1.SetBorderSize(2)
        c1.SetTickx(1)
        c1.SetTicky(1)
        c1.SetLeftMargin(0.15)
        c1.SetRightMargin(0.05)
        c1.SetTopMargin(0.05)
        c1.SetBottomMargin(0.13)
        c1.SetFrameFillStyle(0)
        c1.SetFrameBorderMode(0)
        c1.SetFrameFillStyle(0)
        c1.SetFrameBorderMode(0)
        

        dataPlot.SetFillColor(1)
        dataPlot.SetLineColor(1)
        dataPlot.SetFillStyle(0)
        dataPlot.SetMarkerStyle(20)
        dataPlot.GetXaxis().SetTitle("m_{4#font[12]{l}} (GeV)")
        dataPlot.GetXaxis().SetRange(3,40)
        dataPlot.GetXaxis().SetLabelFont(42)
        dataPlot.GetXaxis().SetLabelSize(0.05)
        dataPlot.GetXaxis().SetTitleSize(0.06)
        dataPlot.GetXaxis().SetTitleOffset(0.9)
        dataPlot.GetXaxis().SetTitleFont(42)
        dataPlot.GetYaxis().SetTitle("Events / 20 GeV")
        dataPlot.GetYaxis().SetLabelFont(42)
        dataPlot.GetYaxis().SetLabelSize(0.05)
        dataPlot.GetYaxis().SetTitleSize(0.06)
        dataPlot.GetYaxis().SetTitleOffset(1.25)
        dataPlot.GetYaxis().SetTitleFont(42)
        dataPlot.GetZaxis().SetLabelFont(42)
        dataPlot.GetZaxis().SetLabelSize(0.035)
        dataPlot.GetZaxis().SetTitleSize(0.035)
        dataPlot.GetZaxis().SetTitleFont(42)
        
        dataPlot.Draw("e1 goff")
        dataPlot.GetXaxis().SetRangeUser(50,800)
        dataPlot.GetYaxis().SetRangeUser(0, 1.3 * dataPlot.GetMaximum())
        A.Draw("hist same goff")
        dataPlot.Draw("same e1 goff")
        if (p > 4):
        
            _fileD = ROOT.TFile("estimateZXData.root")
            histPlot1 = _fileD.Get("h1D_m4l_Add_"+dataName[p-5])
            print ("h1D_m4l_Add_"+dataName[p-5])
            print (histPlot1.Integral())
            leg.AddEntry(histPlot1,"2P2F contr smooth","F")
            histPlot1.SetFillColor(0)
            histPlot1.SetLineColor(6)
            histPlot1.SetLineWidth(2)
            histPlot1.SetFillStyle(4000)
            histPlot1.Smooth()
            histPlot1.Draw("hist same goff")
            

        leg.Draw("goff")
        
        tex1 = ROOT.TLatex(0.95,0.96,"58.8 fb^{-1} (13 TeV)")
        tex1.SetNDC()
        tex1.SetTextAlign(31)
        tex1.SetTextFont(42)
        tex1.SetTextSize(0.03)
        tex1.SetLineWidth(2)
        tex1.Draw("goff")
        tex2 = ROOT.TLatex(0.15,0.96,"CMS")
        tex2.SetNDC()
        tex2.SetTextFont(61)
        tex2.SetTextSize(0.0375)
        tex2.SetLineWidth(2)
        tex2.Draw("goff")
        tex3 = ROOT.TLatex(0.23,0.96,"Preliminary")
        tex3.SetNDC()
        tex3.SetTextFont(52)
        tex3.SetTextSize(0.0285)
        tex3.SetLineWidth(2)
        tex3.Draw("goff")
        
        c1.SaveAs("CR_plots/"+histName[p]+"_match_auto.pdf")
        c1.SaveAs("CR_plots/"+histName[p]+"_match_auto.C")
        c1.SaveAs("CR_plots/"+histName[p]+"_match_auto.root")


def PlotNDProcess(Type, VarP):

    histName = ["el_EB","el_EE","mu_EB","mu_EE"]
    histNameLabel = ["electron barrel","electron endcap","muon barrel","muon endcap"]
    Fillcolors = ["#cc0099","#99ccff","#996666","#669966"]
    Linecolors = ["#990066","#000099","#5f3f3f","#003300"]
    StackTitle = "Proba"
    ProcessNames = ["Data","WZ","ZZ","TT","DY50"]
    fileSize = 5
    plotSize = 4
    variableN = "ptl3"

    for  p in range (plotSize): 
        
        leg_xl = 0.58
        leg_xr = 0.98
        leg_yb = 0.55
        leg_yt = 0.85

        if (VarP == "mZ1"):
            leg_xl = 0.25
            leg_xr = 0.75
            leg_yb = 0.55
            leg_yt = 0.85 



        h = ""

        if (Type == "n"):
            h ="numerator"
        else:
            h = "denominator"
        
        leg = ROOT.TLegend(leg_xl,leg_yb,leg_xr,leg_yt)

        leg.SetBorderSize(0)
        leg.SetLineColor(1)
        leg.SetLineStyle(1)
        leg.SetLineWidth(1)
        leg.SetFillColor(0)
        leg.SetFillStyle(0)

        entry = leg.AddEntry(ROOT.nullptr, histNameLabel[p]+" "+h,"h")
        
        entry.SetLineColor(1)
        entry.SetLineStyle(1)
        entry.SetLineWidth(1)
        entry.SetMarkerColor(1)
        entry.SetMarkerStyle(21)
        entry.SetMarkerSize(1)
        entry.SetTextFont(62)
        entry=leg.AddEntry(ROOT.nullptr,"Data","LP")
        entry.SetLineColor(1)
        entry.SetLineStyle(1)
        entry.SetLineWidth(1)
        entry.SetMarkerColor(1)
        entry.SetMarkerStyle(20)
        entry.SetMarkerSize(1)
        entry.SetTextFont(42)
        entry=leg.AddEntry(ROOT.nullptr,"Z + jets","F")
    
        ci = ROOT.TColor.GetColor("#669966")
        entry.SetFillColor(ci)
        entry.SetFillStyle(1001)
        
        ci = ROOT.TColor.GetColor("#003300")
        entry.SetLineColor(ci)
        entry.SetLineStyle(1)
        entry.SetLineWidth(1)
        entry.SetMarkerColor(1)
        entry.SetMarkerStyle(21)
        entry.SetMarkerSize(1)
        entry.SetTextFont(42)
        entry=leg.AddEntry(ROOT.nullptr,"t#bar{t} + jets","F")
        
        ci = ROOT.TColor.GetColor("#996666")
        entry.SetFillColor(ci)
        entry.SetFillStyle(1001)
        
        ci = ROOT.TColor.GetColor("#5f3f3f")
        entry.SetLineColor(ci)
        entry.SetLineStyle(1)
        entry.SetLineWidth(1)
        entry.SetMarkerColor(1)
        entry.SetMarkerStyle(21)
        entry.SetMarkerSize(1)
        entry.SetTextFont(42)

        entry=leg.AddEntry(ROOT.nullptr,"WZ","F")
        
        ci = ROOT.TColor.GetColor("#cc0099")
        entry.SetFillColor(ci)
        entry.SetFillStyle(1001)
        
        ci = ROOT.TColor.GetColor("#990066")
        entry.SetLineColor(ci)
        entry.SetLineStyle(1)
        entry.SetLineWidth(1)
        entry.SetMarkerColor(1)
        entry.SetMarkerStyle(21)
        entry.SetMarkerSize(1)
        entry.SetTextFont(42)

        entry=leg.AddEntry(ROOT.nullptr,"Z#gamma*,ZZ","F")
        
        ci = ROOT.TColor.GetColor("#99ccff")
        entry.SetFillColor(ci)
        entry.SetFillStyle(1001)
    
        
        A = ROOT.THStack(StackTitle,StackTitle)
        _file = []
        for i in range (fileSize):
            
            Title = "MC"
            if(i == 0):
                Title = "Data"
            
            _file.append(ROOT.TFile("Hist_"+Title+"_"+VarP+"_"+ProcessNames[i]+".root"))
            print (_file)
            if (i==0):
                dataPlot = _file[i].Get("Data_FR"+histName[p]+"_"+Type)
                dataPlot.SetFillColor(1)
                dataPlot.SetLineColor(1)
                dataPlot.SetFillStyle(0)
                dataPlot.SetMarkerStyle(20)
                dataPlot.GetXaxis().SetRange(3,40)
                dataPlot.GetXaxis().SetLabelFont(42)
                dataPlot.GetXaxis().SetLabelSize(0.05)
                dataPlot.GetXaxis().SetTitleSize(0.06)
                dataPlot.GetXaxis().SetTitleOffset(0.9)
                dataPlot.GetXaxis().SetTitleFont(42)
                dataPlot.GetYaxis().SetLabelFont(42)
                dataPlot.GetYaxis().SetLabelSize(0.05)
                dataPlot.GetYaxis().SetTitleSize(0.06)
                dataPlot.GetYaxis().SetTitleOffset(1.6)
                dataPlot.GetYaxis().SetTitleFont(42)
                dataPlot.GetZaxis().SetLabelFont(40)
                dataPlot.GetZaxis().SetLabelSize(0.035)
                dataPlot.GetZaxis().SetTitleSize(0.035)
                dataPlot.GetZaxis().SetTitleFont(42)
                dataPlot.GetXaxis().SetTitle(VarP+" [GeV]")
                
                if (VarP=="mZ1"):
                    dataPlot.GetYaxis().SetTitle("Events/3 GeV")
                if (VarP=="mZ1"):
                    dataPlot.GetXaxis().SetTitle("m_{Z_{1}} (GeV)")
                if (VarP=="mEt"):
                    dataPlot.GetYaxis().SetTitle("Events/5 GeV")
                if (VarP=="mEt"):
                    dataPlot.GetXaxis().SetTitle("E_{T_{miss}} (GeV)")
                if (VarP=="ptl3"):
                    dataPlot.GetXaxis().SetRangeUser(5,80)
                    dataPlot.GetYaxis().SetRangeUser(0,1.1*dataPlot.GetMaximum())                
                if (VarP=="mEt"):
                    dataPlot.GetXaxis().SetRangeUser(0,50)
                    dataPlot.GetYaxis().SetRangeUser(0,1.5*dataPlot.GetMaximum())
                
                if (VarP=="mZ1"):
                    dataPlot.GetXaxis().SetRangeUser(60,120)
                    dataPlot.GetYaxis().SetRangeUser(0,1.5*dataPlot.GetMaximum())
 
            else:
                histPlot = _file[i].Get("Data_FR"+histName[p]+"_"+Type)
                histPlot.SetFillColor(ROOT.TColor.GetColor(Fillcolors[i-1]))
                histPlot.SetLineColor(ROOT.TColor.GetColor(Linecolors[i-1]))
                histPlot.SetFillStyle(1001)
                
                histPlot.SetLineStyle(0)
                histPlot.SetMarkerStyle(20)
                histPlot.GetXaxis().SetLabelFont(42)
                histPlot.GetXaxis().SetLabelOffset(0.007)
                histPlot.GetXaxis().SetLabelSize(0.05)
                histPlot.GetXaxis().SetTitleSize(0.06)
                histPlot.GetXaxis().SetTitleOffset(0.9)
                histPlot.GetXaxis().SetTitleFont(42)
                histPlot.GetYaxis().SetLabelFont(42)
                histPlot.GetYaxis().SetLabelOffset(0.007)
                histPlot.GetYaxis().SetLabelSize(0.05)
                histPlot.GetYaxis().SetTitleSize(0.06)
                histPlot.GetYaxis().SetTitleOffset(1.25)
                histPlot.GetYaxis().SetTitleFont(42)
                histPlot.GetZaxis().SetLabelFont(42)
                histPlot.GetZaxis().SetLabelOffset(0.007)
                histPlot.GetZaxis().SetLabelSize(0.05)
                histPlot.GetZaxis().SetTitleSize(0.06)
                histPlot.GetZaxis().SetTitleFont(42)
                A.Add(histPlot)

        
        # setup environment & canvas
        c1 = ROOT.TCanvas("c1", "myPlots",0,75,600,600)
        ROOT.gStyle.SetOptFit(1)
        ROOT.gStyle.SetOptStat(0)
        ROOT.gStyle.SetOptTitle(0)
        c1.Range(-102.5,-10.38415,847.5,69.4939)
        c1.SetFillColor(0)
        c1.SetBorderMode(0)
        c1.SetBorderSize(2)
        c1.SetTickx(1)
        c1.SetTicky(1)
        c1.SetLeftMargin(0.18)
        c1.SetRightMargin(0.05)
        c1.SetTopMargin(0.05)
        c1.SetBottomMargin(0.13)
        c1.SetFrameFillStyle(0)
        c1.SetFrameBorderMode(0)
        c1.SetFrameFillStyle(0)
        c1.SetFrameBorderMode(0)
        dataPlot.SetMaximum(1.5*dataPlot.GetMaximum())
        if (p>3):
            dataPlot.SetMaximum(2*dataPlot.GetMaximum())

        dataPlot.Draw("e1 goff")

        if (VarP == "mZ1"):
            dataPlot.GetXaxis().SetRangeUser(60,120)
        elif (VarP == "mEt"):
            dataPlot.GetXaxis().SetRangeUser(0,50)

        A.Draw("hist same goff")
        dataPlot.Draw("same e1 goff")
        leg.Draw("goff")
        
        tex = ROOT.TLatex(0.95,0.96,"7.7 fb^{-1} (13 TeV)")
        tex.SetNDC()
        tex.SetTextAlign(31)
        tex.SetTextFont(42)
        tex.SetTextSize(0.03)
        tex.SetLineWidth(2)
        tex.Draw("goff")
        tex1 = ROOT.TLatex(0.15,0.96,"CMS")
        tex1.SetNDC()
        tex1.SetTextFont(61)
        tex1.SetTextSize(0.0375)
        tex1.SetLineWidth(2)
        tex1.Draw("goff")
        tex2 = ROOT.TLatex(0.23,0.96,"Preliminary")
        tex2.SetNDC()
        tex2.SetTextFont(52)
        tex2.SetTextSize(0.0285)
        tex2.SetLineWidth(2)
        tex2.Draw("goff")
        
        
        c1.SaveAs("Z1L_plots/"+VarP+histName[p]+"_"+Type+"_process_auto.pdf")

def PlotNDMatched(Type, VarP):
    
    histNameLabel = ["electron endcap","electron barrel","muon endcap","muon barrel"]
    
    histName = ["EE","EB","ME","MB"]
    dataName = ["el_EE","el_EB","mu_EE","mu_EB"]
    StackTitle = "Proba"
    
    fileSize = 5
    plotSize = 4
    variableN = "ptl3"

    for p in range (plotSize):
        
        h = ""
        if (Type == "n"):
            h ="numerator"
        else:
            h = "denominator"
        

        leg_xl = 0.58
        leg_xr = 0.98
        leg_yb = 0.55
        leg_yt = 0.85

        if (VarP == "mZ1"):
            leg_xl = 0.25
            leg_xr = 0.75
            leg_yb = 0.55
            leg_yt = 0.85 

        leg = ROOT.TLegend(leg_xl,leg_yb,leg_xr,leg_yt)
        
        _file1 = ROOT.TFile("Hist_Data_"+VarP+"_Data.root")

        print ("Data_FR"+dataName[p]+"_"+Type)
        
        dataPlot = _file1.Get("Data_FR"+dataName[p]+"_"+Type)
        
        dataPlot.SetFillColor(1)
        dataPlot.SetLineColor(1)
        dataPlot.SetFillStyle(0)
        dataPlot.SetMarkerStyle(20)
        dataPlot.GetXaxis().SetRange(3,40)
        dataPlot.GetXaxis().SetLabelFont(42)
        dataPlot.GetXaxis().SetLabelSize(0.05)
        dataPlot.GetXaxis().SetTitleSize(0.06)
        dataPlot.GetXaxis().SetTitleOffset(0.9)
        dataPlot.GetXaxis().SetTitleFont(42)
        dataPlot.GetYaxis().SetLabelFont(42)
        dataPlot.GetYaxis().SetLabelSize(0.05)
        dataPlot.GetYaxis().SetTitleSize(0.06)
        dataPlot.GetYaxis().SetTitleOffset(1.6)
        dataPlot.GetYaxis().SetTitleFont(42)
        dataPlot.GetZaxis().SetLabelFont(40)
        dataPlot.GetZaxis().SetLabelSize(0.035)
        dataPlot.GetZaxis().SetTitleSize(0.035)
        dataPlot.GetZaxis().SetTitleFont(42)
        dataPlot.GetXaxis().SetTitle(VarP+" [GeV]")
        
        if (VarP=="mZ1"):
            dataPlot.GetYaxis().SetTitle("Events/3 GeV")
        if (VarP=="mZ1"):
            dataPlot.GetXaxis().SetTitle("m_{Z_{1}} (GeV)")
        if (VarP=="mEt"):
            dataPlot.GetYaxis().SetTitle("Events/5 GeV")
        if (VarP=="mEt"):
            dataPlot.GetXaxis().SetTitle("E_{T_{miss}} (GeV)")
        if (VarP=="ptl3"):
            dataPlot.GetXaxis().SetRangeUser(5,80)
            dataPlot.GetYaxis().SetRangeUser(0,1.1*dataPlot.GetMaximum())                
        if (VarP=="mEt"):
            dataPlot.GetXaxis().SetRangeUser(0,50)
            dataPlot.GetYaxis().SetRangeUser(0,1.5*dataPlot.GetMaximum())
        
        if (VarP=="mZ1"):
            dataPlot.GetXaxis().SetRangeUser(60,120)
            dataPlot.GetYaxis().SetRangeUser(0,1.5*dataPlot.GetMaximum())


        if (VarP=="mZ1"):
            dataPlot.GetYaxis().SetTitle("Events/3 GeV")
            dataPlot.GetXaxis().SetTitle("m_{Z_{1}} (GeV)")
        if (VarP=="mEt"):
            dataPlot.GetYaxis().SetTitle("Events/5 GeV")
            dataPlot.GetXaxis().SetTitle("E_{T_{miss}} (GeV)")
        

        Title = "Data"
        leg.AddEntry(dataPlot,Title,"LP")
        
        _fileMC = ROOT.TFile("Hist_ID_"+VarP+".root")
        A = ROOT.THStack(StackTitle,StackTitle)
        
        histPlotC = _fileMC.Get("Hist_conv"+histName[p]+"_"+Type)
        histPlotC.SetFillColor(30)
        histPlotC.SetLineColor(1)
        histPlotC.SetFillStyle(1001)
        
        histPlotC.SetLineStyle(0)
        histPlotC.SetMarkerStyle(20)
        histPlotC.GetXaxis().SetLabelFont(42)
        histPlotC.GetXaxis().SetLabelOffset(0.007)
        histPlotC.GetXaxis().SetLabelSize(0.05)
        histPlotC.GetXaxis().SetTitleSize(0.06)
        histPlotC.GetXaxis().SetTitleOffset(0.9)
        histPlotC.GetXaxis().SetTitleFont(42)
        histPlotC.GetYaxis().SetLabelFont(42)
        histPlotC.GetYaxis().SetLabelOffset(0.007)
        histPlotC.GetYaxis().SetLabelSize(0.05)
        histPlotC.GetYaxis().SetTitleSize(0.06)
        histPlotC.GetYaxis().SetTitleOffset(1.25)
        histPlotC.GetYaxis().SetTitleFont(42)
        histPlotC.GetZaxis().SetLabelFont(42)
        histPlotC.GetZaxis().SetLabelOffset(0.007)
        histPlotC.GetZaxis().SetLabelSize(0.05)
        histPlotC.GetZaxis().SetTitleSize(0.06)
        histPlotC.GetZaxis().SetTitleFont(42)
        
        
        leg.AddEntry(histPlotC,"Conversions","F")
        A.Add(histPlotC)
        
        histPlotP = _fileMC.Get("Hist_prompt"+histName[p]+"_"+Type)
        histPlotP.SetFillColor(38)
        histPlotP.SetLineColor(1)
        histPlotP.SetFillStyle(1001)
        histPlotP.SetLineStyle(0)
        histPlotP.SetMarkerStyle(20)
        histPlotP.GetXaxis().SetLabelFont(42)
        histPlotP.GetXaxis().SetLabelOffset(0.007)
        histPlotP.GetXaxis().SetLabelSize(0.05)
        histPlotP.GetXaxis().SetTitleSize(0.06)
        histPlotP.GetXaxis().SetTitleOffset(0.9)
        histPlotP.GetXaxis().SetTitleFont(42)
        histPlotP.GetYaxis().SetLabelFont(42)
        histPlotP.GetYaxis().SetLabelOffset(0.007)
        histPlotP.GetYaxis().SetLabelSize(0.05)
        histPlotP.GetYaxis().SetTitleSize(0.06)
        histPlotP.GetYaxis().SetTitleOffset(1.25)
        histPlotP.GetYaxis().SetTitleFont(42)
        histPlotP.GetZaxis().SetLabelFont(42)
        histPlotP.GetZaxis().SetLabelOffset(0.007)
        histPlotP.GetZaxis().SetLabelSize(0.05)
        histPlotP.GetZaxis().SetTitleSize(0.06)
        histPlotP.GetZaxis().SetTitleFont(42)
        
        leg.AddEntry(histPlotP,"Prompt","F")
        A.Add(histPlotP)
        
        histPlotF = _fileMC.Get("Hist_fakes"+histName[p]+"_"+Type)
        histPlotF.SetFillColor(46)
        histPlotF.SetLineColor(1)
        histPlotF.SetFillStyle(1001)
        histPlotF.SetLineStyle(0)
        histPlotF.SetMarkerStyle(20)
        histPlotF.GetXaxis().SetLabelFont(42)
        histPlotF.GetXaxis().SetLabelOffset(0.007)
        histPlotF.GetXaxis().SetLabelSize(0.05)
        histPlotF.GetXaxis().SetTitleSize(0.06)
        histPlotF.GetXaxis().SetTitleOffset(0.9)
        histPlotF.GetXaxis().SetTitleFont(42)
        histPlotF.GetYaxis().SetLabelFont(42)
        histPlotF.GetYaxis().SetLabelOffset(0.007)
        histPlotF.GetYaxis().SetLabelSize(0.05)
        histPlotF.GetYaxis().SetTitleSize(0.06)
        histPlotF.GetYaxis().SetTitleOffset(1.25)
        histPlotF.GetYaxis().SetTitleFont(42)
        histPlotF.GetZaxis().SetLabelFont(42)
        histPlotF.GetZaxis().SetLabelOffset(0.007)
        histPlotF.GetZaxis().SetLabelSize(0.05)
        histPlotF.GetZaxis().SetTitleSize(0.06)
        histPlotF.GetZaxis().SetTitleFont(42)
        
        leg.AddEntry(histPlotF,"LF Fakes","F")
        A.Add(histPlotF)
        
        histPlotB = _fileMC.Get("Hist_BDfakes"+histName[p]+"_"+Type)
        histPlotB.SetFillColor(49)
        histPlotB.SetLineColor(1)
        histPlotB.SetFillStyle(1001)
        histPlotB.SetLineStyle(0)
        histPlotB.SetMarkerStyle(20)
        histPlotB.GetXaxis().SetLabelFont(42)
        histPlotB.GetXaxis().SetLabelOffset(0.007)
        histPlotB.GetXaxis().SetLabelSize(0.05)
        histPlotB.GetXaxis().SetTitleSize(0.06)
        histPlotB.GetXaxis().SetTitleOffset(0.9)
        histPlotB.GetXaxis().SetTitleFont(42)
        histPlotB.GetYaxis().SetLabelFont(42)
        histPlotB.GetYaxis().SetLabelOffset(0.007)
        histPlotB.GetYaxis().SetLabelSize(0.05)
        histPlotB.GetYaxis().SetTitleSize(0.06)
        histPlotB.GetYaxis().SetTitleOffset(1.25)
        histPlotB.GetYaxis().SetTitleFont(42)
        histPlotB.GetZaxis().SetLabelFont(42)
        histPlotB.GetZaxis().SetLabelOffset(0.007)
        histPlotB.GetZaxis().SetLabelSize(0.05)
        histPlotB.GetZaxis().SetTitleSize(0.06)
        histPlotB.GetZaxis().SetTitleFont(42)
        
        leg.AddEntry(histPlotB,"HF fakes","F")
        A.Add(histPlotB)
        
        # setup environment & canvas
        c1 = ROOT.TCanvas("c1", "myPlots",0,75,600,600)
        ROOT.gStyle.SetOptFit(1)
        ROOT.gStyle.SetOptStat(0)
        ROOT.gStyle.SetOptTitle(0)
        c1.Range(-102.5,-10.38415,847.5,69.4939)
        c1.SetFillColor(0)
        c1.SetBorderMode(0)
        c1.SetBorderSize(2)
        c1.SetTickx(1)
        c1.SetTicky(1)
        c1.SetLeftMargin(0.18)
        c1.SetRightMargin(0.05)
        c1.SetTopMargin(0.05)
        c1.SetBottomMargin(0.13)
        c1.SetFrameFillStyle(0)
        c1.SetFrameBorderMode(0)
        c1.SetFrameFillStyle(0)
        c1.SetFrameBorderMode(0)
        dataPlot.Draw("e1")

        dataPlot.SetMaximum(1.5*dataPlot.GetMaximum())
        if (p>3):
            dataPlot.SetMaximum(2*dataPlot.GetMaximum())

        if (VarP == "mZ1"):
            dataPlot.GetXaxis().SetRangeUser(60,120)
        elif (VarP == "mEt"):
            dataPlot.GetXaxis().SetRangeUser(0,50)

        dataPlot.Draw("e1 goff")

        A.Draw("hist same goff")
        dataPlot.Draw("same e1 goff")
        leg.Draw("goff")
        
        
        tex = ROOT.TLatex(0.95,0.96,"58.8 fb^{-1} (13 TeV)")
        tex.SetNDC()
        tex.SetTextAlign(31)
        tex.SetTextFont(42)
        tex.SetTextSize(0.03)
        tex.SetLineWidth(2)
        tex.Draw("goff")
        tex1 = ROOT.TLatex(0.15,0.96,"CMS")
        tex1.SetNDC()
        tex1.SetTextFont(61)
        tex1.SetTextSize(0.0375)
        tex1.SetLineWidth(2)
        tex1.Draw("goff")
        tex2 = ROOT.TLatex(0.23,0.96,"Preliminary")
        tex2.SetNDC()
        tex2.SetTextFont(52)
        tex2.SetTextSize(0.0285)
        tex2.SetLineWidth(2)
        tex2.Draw("goff")
        
        
        c1.SaveAs("Z1L_plots/"+VarP+"_"+Type+"_"+histName[p]+"_match_auto.pdf")
    


os.system("mkdir CR_plots")
os.system("mkdir Z1L_plots")
#CRPlot()
#CRMatchPlot()
#PlotNDProcess("n","ptl3")
#PlotNDProcess("d","ptl3")
PlotNDProcess("n","mZ1")
PlotNDProcess("d","mZ1")
PlotNDProcess("n","mEt")
PlotNDProcess("d","mEt")
#PlotNDMatched("n","ptl3")
#PlotNDMatched("d","ptl3")