from ROOT import *
import sys,os
from ROOT import TStopwatch
c1 = TCanvas()
c1.SetGridx()
c1.SetGridy()
c1.SetTickx()
c1.SetTicky()

#_________________________________________________________________________________________

def MidPointIntegral(a = 0, b = 1, n=1000, func = TF1("f1", " x ", 0, 1)):
    deltax = (1.0*(1.0*b-1.0*a)) / (1.0*n)
    int1 = 0.0
    for x in range (0,n):
        x1 = 1.0*a + 1.0*deltax*(1.0*x + 0.5) 
        y = 1.0*func.Eval(x1, 0, 0)
        int1 = int1 + y*deltax
    return int1
#_________________________________________________________________________________________

def TrapezoidalIntegral(a = 0, b=1, n = 1000, func = TF1("f1", " x ", 0, 1)):
    deltax = 1.0*(1.0*b-1.0*a) / (1.0*n)
    int1 = 0.0
    for x in range (0.,n):
        x1 = 1.0*a + deltax*x
        x2 = 1.0*a + deltax*(x+1)
        y1 = 1.0*func.Eval(x1, 0, 0)
        y2 = 1.0*func.Eval(x2, 0, 0)
        int1 = int1 + (y1+y2)*0.5*deltax
    return int1
#______________________________________________________________________________________
# Change values below.        

a = 0.                              # Lower Integration Limit
b = 3.                              # Higher Integration Limit

ff = "300*exp(-100000*(x-2)*(x-2)) + 3*exp(x-2)"
func = TF1("f1", ff, a, b)

ww = "1"
ww = "exp(-1*(x-2)*(x-2)) "
ww = "exp(x-2) + exp(-10000*(x-2)*(x-2)) "

wht = TF1("ww", ww, a,b)

minnumberofsteps = 1                # Min number of divisions
maxnumberofsteps = 1000             # Max number of divisions

trsh = 0.001                        # Point of convergence for MC
probe_mc_every = 100

#_________________________________________________________________________________________

timer_Midpoint = TH1F("timer_Midpoint","timer_Midpoint",maxnumberofsteps, minnumberofsteps,maxnumberofsteps)
timer_Midpoint.SetLineColor(kBlue)
timer_Midpoint.SetLineWidth(4)
Integ_Midpoint = TH1F("Integ_Midpoint","Integ_Midpoint",maxnumberofsteps, minnumberofsteps,maxnumberofsteps)
Integ_Midpoint.SetLineColor(kBlue)
Integ_Midpoint.SetLineWidth(6)
Integ_Err_Midpoint= TH1F("Integ_Err_Midpoint","Integ_Err_Midpoint",maxnumberofsteps, minnumberofsteps,maxnumberofsteps)
Integ_Err_Midpoint.SetLineColor(kBlue)
Integ_Err_Midpoint.SetLineWidth(3)

timer_1 = TStopwatch()
timer_Cpu_Midpoint = 0.0
Midpoint_stops_at_n = 0.0
old_integration_value = 0.0
first_time_flag = 0

for x in range (minnumberofsteps, maxnumberofsteps+1):
   timer_1.Start()
   Integ_Midpoint.SetBinContent(x, MidPointIntegral(a, b, x, func))             # Set the value of the integral as a function of the steps x
   y = old_integration_value - 1.0*MidPointIntegral(a, b, x, func)              # Calculate Delta I
   y = 1.0*abs(y)/MidPointIntegral(a, b, x, func)                               # Calculate Delta I/I
   Integ_Err_Midpoint.SetBinContent(x, y)                                       # Set the Error histogram equal to Delta I / I
   old_integration_value = MidPointIntegral(a, b, x, func)                      # Store the "old" integration value
   timer_1.Stop()
                
   timer_Cpu_Midpoint = timer_Cpu_Midpoint + timer_1.CpuTime()
   #timer_Real_Midpoint = timer_Real_Midpoint + timer_1.RealTime()
   timer_Midpoint.SetBinContent(x, timer_Cpu_Midpoint) 
   
   if y< trsh and first_time_flag == 0:
        first_time_flag = 1
        Midpoint_stops_at_n = x
        print "Midpoint Done at: ",x
        print " Integral value MidPoint: ",MidPointIntegral(a, b, x, func)
        print " Error value MidPoint: ",y
     
timer_2 = TStopwatch()
timer_Cpu_Trapezoid = 0.0
Trapezoid_stops_at_n = 0.0
old_integration_value = 0.0
first_time_flag = 0

timer_Trapezoid = TH1F("timer_Trapezoid ", "timer_Trapezoid ", maxnumberofsteps, minnumberofsteps,maxnumberofsteps)
timer_Trapezoid.SetLineColor(kRed)
timer_Trapezoid.SetLineWidth(4)
Integ_Trapezoid = TH1F("Integ_Trapezoid", "Integ_Trapezoid", maxnumberofsteps, minnumberofsteps,maxnumberofsteps)
Integ_Trapezoid.SetLineColor(kRed)
Integ_Trapezoid.SetLineWidth(3)
Integ_Err_Trapezoid = TH1F("Integ_Err_Trapezoid", "Integ_Err_Trapezoid", maxnumberofsteps, minnumberofsteps,maxnumberofsteps)
Integ_Err_Trapezoid.SetLineColor(kRed)
Integ_Err_Trapezoid.SetLineWidth(3)

for x in range (minnumberofsteps, maxnumberofsteps+1):
    timer_2.Start()
    Integ_Trapezoid.SetBinContent(x, TrapezoidalIntegral(a, b, x, func))
    y = 1.0*TrapezoidalIntegral(a, b, x, func) - old_integration_value
    y = 1.0*abs(y)/TrapezoidalIntegral(a, b, x, func)
    Integ_Err_Trapezoid.SetBinContent(x,y)
    timer_2.Stop()
    old_integration_value = TrapezoidalIntegral(a, b, x, func)
                               
    timer_Cpu_Trapezoid = timer_Cpu_Trapezoid + timer_2.CpuTime()
    #timeracr = timeracr + xac.RealTime()
    timer_Trapezoid.SetBinContent(x, timer_Cpu_Trapezoid)

    if first_time_flag == 0 and abs(y)< trsh:
        first_time_flag = 1
        Trapezoid_stops_at_n = x
        
        print "Trapezoid Done at: ",Trapezoid_stops_at_n
        print " Integral Value Trapezoid: ",TrapezoidalIntegral(a, b, x, func)
        print " Error Value of Trapezoid: ",y
        
         
Integ_MonteCarlo = TH1D("Integ_MonteCarlo ", "Integ_MonteCarlo ", maxnumberofsteps, minnumberofsteps,maxnumberofsteps)
Integ_MonteCarlo.SetLineColor(kBlack)
Integ_MonteCarlo.SetLineWidth(10)
Integ_Err_MonteCarlo = TH1D("Integ_Err_MonteCarlo", "Integ_Err_MonteCarlo", maxnumberofsteps, minnumberofsteps,maxnumberofsteps)
Integ_Err_MonteCarlo.SetLineColor(kBlack)
Integ_Err_MonteCarlo.SetLineWidth(9)
Integ_RMS_MonteCarlo = TH1D("Integ_RMS_MonteCarlo", "Integ_RMS_MonteCarlo", maxnumberofsteps, minnumberofsteps,maxnumberofsteps)
Integ_RMS_MonteCarlo.SetLineColor(kGreen+2)
Integ_RMS_MonteCarlo.SetLineWidth(9)
timer_MonteCarlo = TH1F("mctimer", "mctimer", maxnumberofsteps, minnumberofsteps, maxnumberofsteps)
timer_MonteCarlo.SetLineColor(kBlack)
timer_MonteCarlo.SetLineWidth(5)

timer_3 = TStopwatch()
timer_Cpu_MonteCarlo = 0.0
MonteCarlo_stops_at_n = 0.0
first_time_flag = 0.0

sumf = 0.0
sumf2 = 0.0

for x1 in range (minnumberofsteps, maxnumberofsteps+1):
    timer_3.Start()
    for c in range (0., probe_mc_every+1):                                                                      
        xa = 1.0*a + (1.0*b - 1.0*a)*TRandom1().Rndm()
        sumf = sumf + 1.0*func.Eval(xa,0,0)
        sumf2 = sumf2 + 1.0*func.Eval(xa,0,0)*1.0*func.Eval(xa,0,0)                       

    mcintegral = (1.0*b - 1.0*a)*sumf/probe_mc_every   
    mc_rms = abs(sumf2 - sumf*sumf)/((x1*probe_mc_every)*(x1*probe_mc_every))    
    Integ_MonteCarlo.SetBinContent(int(x1), mcintegral/x1)
    Integ_RMS_MonteCarlo.SetBinContent(x1, sqrt(mc_rms))
    Integ_Err_MonteCarlo.SetBinContent(x1, sqrt(mc_rms)/(x1*(probe_mc_every)))
    timer_3.Stop()
    
    timer_Cpu_MonteCarlo= timer_Cpu_MonteCarlo + timer_3.CpuTime()
    timer_MonteCarlo.SetBinContent(x1, timer_Cpu_MonteCarlo)

    if first_time_flag == 0 and abs(sqrt(mc_rms)/(x1*(probe_mc_every)))< trsh:
        first_time_flag = 1
        MonteCarlo_stops_at_n = x1
        
        print "MCDone at: ",MonteCarlo_stops_at_n
        print " Integral Value MC: ",mcintegral/x1
        print " Error Value of MC: ",abs(sqrt(mc_rms)/(x1*(probe_mc_every)))
        
        
Integ_MonteCarlo_Weight = TH1D("Integ_MonteCarlo_Weight ", "Integ_MonteCarlo_Weight ", maxnumberofsteps, minnumberofsteps,maxnumberofsteps)
Integ_MonteCarlo_Weight.SetLineColor(kOrange)
Integ_MonteCarlo_Weight.SetLineWidth(4)
Integ_Err_MonteCarlo_Weight = TH1D("Integ_Err_MonteCarlo_Weight", "Integ_Err_MonteCarlo_Weight", maxnumberofsteps, minnumberofsteps,maxnumberofsteps)
Integ_Err_MonteCarlo_Weight.SetLineColor(kOrange)
Integ_Err_MonteCarlo_Weight.SetLineWidth(4)
Integ_RMS_MonteCarlo_Weight = TH1D("Integ_RMS_MonteCarlo_Weight", "Integ_RMS_MonteCarlo_Weight", maxnumberofsteps, minnumberofsteps,maxnumberofsteps)
Integ_RMS_MonteCarlo_Weight.SetLineColor(kMagenta)
Integ_RMS_MonteCarlo_Weight.SetLineWidth(4)
timer_MonteCarlo_Weight = TH1F("mctimer_Weight", "mctimer_Weight", maxnumberofsteps, minnumberofsteps, maxnumberofsteps)
timer_MonteCarlo_Weight.SetLineColor(kOrange)
timer_MonteCarlo_Weight.SetLineWidth(6)


timer_4 = TStopwatch()
timer_Cpu_MonteCarlo_Weighted  = 0.0
MonteCarlo_Weighted_stops_at_n = 0.0
first_time_flag = 0.0

sumf = 0.0
sumf2 = 0.0
weh = 0.0

for x1 in range (minnumberofsteps, maxnumberofsteps+1):
    timer_4.Start()
    for c in range (0., probe_mc_every+1):

        xa = wht.GetRandom(a,b)
        weh = weh + 1.0/(1.0*wht.Eval(xa, 0, 0))
        sumf = sumf + 1.0*func.Eval(xa,0,0)/(1.0*wht.Eval(xa, 0, 0))
        sumf2 = sumf2 + 1.0*func.Eval(xa,0,0)*1.0*func.Eval(xa,0,0) /(1.0*wht.Eval(xa, 0, 0)*wht.Eval(xa, 0, 0))                      



    mcintegral = (1.0*b - 1.0*a)*sumf/weh   
    mc_rms = abs(sumf2 - sumf*sumf)/((weh)*(weh))

    
    Integ_MonteCarlo_Weight.SetBinContent(int(x1), mcintegral)
    Integ_RMS_MonteCarlo_Weight.SetBinContent(x1, sqrt(mc_rms))
    Integ_Err_MonteCarlo_Weight.SetBinContent(x1, sqrt(mc_rms)/((weh)))


    timer_4.Stop()
    timer_Cpu_MonteCarlo_Weighted= timer_Cpu_MonteCarlo_Weighted + timer_4.CpuTime()
    timer_MonteCarlo_Weight.SetBinContent(x1, timer_Cpu_MonteCarlo_Weighted)




    if first_time_flag == 0 and abs(sqrt(mc_rms)/(x1*(probe_mc_every)))< trsh:
        first_time_flag = 1
        MonteCarlo_Weighted_stops_at_n = x1
        
        print "MC Weighted Done at: ",MonteCarlo_Weighted_stops_at_n
        print " Integral Value MC: ",mcintegral
        print " Error Value of MC: ",abs(sqrt(mc_rms)/(x1*(probe_mc_every)))








maxy = Integ_Trapezoid.GetMaximum()
miny = 0

b1 = TH2F("b1","b1",maxnumberofsteps,0,maxnumberofsteps,50,miny, 0.1*maxy)
text = TLatex(maxnumberofsteps*0.2, maxy*0.75, " f(x) = "+ff)
text2 = TLatex(maxnumberofsteps*0.2, maxy*0.62, " %s < x < %s"%(a,b))



########################## Integral Convergence Comparison ########################333


legend = TLegend(0.52,0.6,0.89,0.84)
legend.SetFillColor(ROOT.kWhite)
legend.SetMargin(0.25)
legend.SetHeader(" Integral Method:")
legend.AddEntry(Integ_Midpoint,"Midpoint", "l")
legend.AddEntry(Integ_Trapezoid,"Trapezoid","l")
legend.AddEntry(Integ_MonteCarlo,"Monte Carlo ","l")
legend.AddEntry(Integ_MonteCarlo_Weight,"Monte Carlo W.","l")



b1.GetXaxis().SetTitle("Number of Steps")
b1.GetYaxis().SetTitle("Integral Value")
b1.SetTitle("Numerical Integration: Convergence Comparison")
b1.SetStats(0)





b1.Draw()
Integ_Trapezoid.Draw("same")
Integ_Midpoint.Draw("same")

Integ_MonteCarlo.Draw("same")
Integ_MonteCarlo_Weight.Draw("same")
legend.Draw("same")

text.Draw("same")
text2.Draw("same")

c1.SaveAs("Lab1_IntegralComparison_W.pdf")



##########################################################################################
maxy = Integ_RMS_MonteCarlo.GetMaximum()*1.5
b1 = TH2F("b1","b1",maxnumberofsteps,0,maxnumberofsteps,50,-0.5, maxy )
b1.GetXaxis().SetTitle("Number of Steps")
b1.GetYaxis().SetTitle("Integral Value")
b1.SetTitle("Monte Carlo  Comparison")
b1.SetStats(0)



legend = TLegend(0.52,0.11,0.89,0.35)
legend.SetFillColor(ROOT.kWhite)
legend.SetMargin(0.25)
legend.SetHeader("Monte Carlo Value:")
legend.AddEntry(Integ_RMS_MonteCarlo,"RMS Value", "l")
legend.AddEntry(Integ_Err_MonteCarlo,"| \Delta I | / I","l")
legend.AddEntry(Integ_RMS_MonteCarlo_Weight,"RMS Value (W)", "l")
legend.AddEntry(Integ_Err_MonteCarlo_Weight,"| \Delta I | / I (W)","l") 

b1.Draw()
Integ_RMS_MonteCarlo.Draw("same")
Integ_Err_MonteCarlo.Draw("same")
Integ_RMS_MonteCarlo_Weight.Draw("same")
Integ_Err_MonteCarlo_Weight.Draw("same")
legend.Draw("same")
text.Draw("same")
text2.Draw("same")
c1.SaveAs("Lab1_MonteCarlo_W.pdf")





############################################################################################

maxy = timer_MonteCarlo.GetMaximum()*2
miny = -0.1
b1 = TH2F("b1","b1",maxnumberofsteps,0,maxnumberofsteps,50,miny, maxy)
b1.GetXaxis().SetTitle("Number of Steps")
b1.GetYaxis().SetTitle("Integral Value")
b1.SetTitle("Numerical Integration: Timer Comparison")
b1.SetStats(0)

legend = TLegend(0.52,0.11,0.89,0.35)
legend.SetFillColor(ROOT.kWhite)
legend.SetMargin(0.25)
legend.SetHeader("Integral Method: ")
legend.AddEntry(timer_MonteCarlo,"Monte Carlo (%d)"%MonteCarlo_stops_at_n, "l")
legend.AddEntry(timer_Trapezoid," Trapezoid (%d)"%Trapezoid_stops_at_n,"l") 
legend.AddEntry(timer_Midpoint," Midpoint (%d)"%Midpoint_stops_at_n,"l") 
legend.AddEntry(timer_MonteCarlo_Weight,"Monte Carlo (%d)"%MonteCarlo_Weighted_stops_at_n, "l")



b1.Draw()
timer_MonteCarlo.Draw("same")
timer_Trapezoid.Draw("same")
timer_Midpoint.Draw("Same")
timer_MonteCarlo_Weight.Draw("same")
legend.Draw("same")
c1.SaveAs("IntegrationTimer_W.pdf")


##################################################################################################
maxy = Integ_Err_Trapezoid.GetMaximum()/5
miny = 1e-9

b1 = TH2F("b1","b1",maxnumberofsteps,0.5,maxnumberofsteps,50,miny, maxy)
b1.GetXaxis().SetTitle("Number of Steps")
b1.GetYaxis().SetTitle("Error ")
b1.SetTitle("Numerical Integration: | \Delta I |/I Comparison")
b1.SetStats(0)

legend = TLegend(0.52,0.11,0.89,0.35)
legend.SetFillColor(ROOT.kWhite)
legend.SetMargin(0.25)
legend.SetHeader("Integral Method: ")
legend.AddEntry(Integ_Err_MonteCarlo,"Monte Carlo (%d)"%MonteCarlo_stops_at_n, "l")
legend.AddEntry(timer_MonteCarlo_Weight,"Monte Carlo (%d)"%MonteCarlo_Weighted_stops_at_n, "l")
legend.AddEntry(Integ_Err_Trapezoid," Trapezoid (%d)"%Trapezoid_stops_at_n,"l") 
legend.AddEntry(Integ_Err_Midpoint," Midpoint (%d)"%Midpoint_stops_at_n,"l")



c1.SetLogy()
#c1.SetLogx()
b1.Draw()

Integ_Err_Trapezoid.Draw("same")
Integ_Err_Midpoint.Draw("same")
Integ_Err_MonteCarlo.Draw("same")
Integ_Err_MonteCarlo_Weight.Draw("same")
legend.Draw("same")
c1.SaveAs("Error_Comparison.pdf")
