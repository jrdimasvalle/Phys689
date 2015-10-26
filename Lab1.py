from ROOT import *
import sys,os
from ROOT import TStopwatch
c1 = TCanvas()
c1.SetGridx()
c1.SetGridy()
c1.SetTickx()
c1.SetTicky()



#_________________________________________________________________________________________

def simpleintegral(a = 0, b = 1, n=1000, func = TF1("f1", " x ", 0, 1)):

    deltax = (1.0*(1.0*b-1.0*a)) / (1.0*n)
    int1 = 0.0
       
    for x in range (0,n):
        x1 = 1.0*a + 1.0*deltax*(1.0*x + 0.5) 
        y = 1.0*func.Eval(x1, 0, 0)
        int1 = int1 + y*deltax


    return int1


#_________________________________________________________________________________________

def trapezoidalint(a = 0, b=1, n = 1000, func = TF1("f1", " x ", 0, 1)):

    deltax = 1.0*(1.0*b-1.0*a) / (1.0*n)
    int1 = 0.0
  
    for x in range (0.,n):
        x1 = 1.0*a + deltax*x
        x2 = 1.0*a + deltax*(x+1)
        y1 = 1.0*func.Eval(x1, 0, 0)
        y2 = 1.0*func.Eval(x2, 0, 0)
        int1 = int1 + (y1+y2)*0.5*deltax
       
    sum2 = 1.0*maxSecondDerivativeValue(a,b,n,func)*(1.0*b - 1.0*a)*(1.0*b - 1.0*a)*(1.0*b - 1.0*a)/(12.0*n*n)


    return int1

#___________________________________________________________________________________________

def MCSimpleIntegral(a = 0, b = 1, n = 1000, func = TF1("f1", "f1", 0 , 1 ), errorx = 0):
        ab = TH1F("ab","ab",n, a,b)
        sum1 = 0.0
        sum2 = 0.0
        
        for c in range (0., n):                                                             # c counts the steps taken from 0 to n
            x1 = 1.0*a + (1.0*b - 1.0*a)*TRandom1().Rndm(1)
            sum1 = sum1 + 1.0*func.Eval(x1,0,0)
            sum2 = sum2 + 1.0*func.Eval(x1,0,0)*1.0*func.Eval(x1,0, 0)                      # Calculate Sum (f(xi) ^2)

        sum1 = sum1 /(1.0*n)                                                                # This is the average of f(x) after the entire set of steps.

        
        sum2 = sum2/(1.0*n) - sum1*sum1                                                     # calculate average of f(xi) ^2
       
        if errorx == 1:
            return abs(1.0*sum2)                                               
        if errorx == 0:
            return sum1*(1.0*b - 1.0*a)                                 
    
#_________________________________________________________________________________________


def MCIntegralWeighted(a = 0, b = 1, n = 1000, func = TF1("f1", " x", 0 , 1 ), wht = TF1("w1", "w1", 0, 1), errorx = 0):
        sum1 = 0.0
        sum2 = 0.0
        whh2 =0.0
        wh = 0.0

        for c in range (0., n+1):                                                          # c counts the steps taken from 0 to n
            x1 = wht.GetRandom(a, b)
            
            sum1 = sum1 + 1.0*func.Eval(x1,0,0)/(1.0*wht.Eval(x1, 0, 0))
            sum2 = sum2 + 1.0*func.Eval(x1,0,0)*1.0*func.Eval(x1,0,0)/(1.0*wht.Eval(x1, 0, 0)*wht.Eval(x1, 0, 0))

                        
            whh2 = whh2 + 1/(1.0*wht.Eval(x1, 0, 0)*wht.Eval(x1, 0, 0))
            wh = wh + 1/(1.0*wht.Eval(x1, 0, 0))


        sum1 = sum1 /(1.0*wh)                       
        sum2 = (sum2/(1.0*whh2))  - sum1*sum1
        
        if errorx ==0:
            return sum1*(1.0*b - 1.0*a)                                 
        if errorx == 1:
            return abs(sum2)


#_________________________________________________________________________________________
# Change values below.        



a = 0.                              # Lower Integration Limit
b = 3.                              # Higher Integration Limit
ww = "300*exp(-100000*(x-2)*(x-2)) + 3*exp(x-2)"
ww = " 1 + sin(x)"
ff = ww
#ff = "1/( (2*x*(x-1) + 1)) "
#ff = "1"
func = TF1("f1", ff, a, b)

#ww = "x/(1 + x) "
weight = TF1("w1", ww, a,b)
minnumberofsteps = 1                # Min number of divisions
maxnumberofsteps = 100             # Max number of divisions

trsh = 0.001                        # Point of convergence for MC

#_________________________________________________________________________________________ 




timerSimple = TH1F("TTsimple","TTsimple",maxnumberofsteps, minnumberofsteps,maxnumberofsteps)
simple = TH1F("simple","simple",maxnumberofsteps, minnumberofsteps,maxnumberofsteps)
simpleerr = TH1F("simpleerr","simpleerr",maxnumberofsteps, minnumberofsteps,maxnumberofsteps)


xtt = TStopwatch()

timerab = 0.0
timerabr = 0.0
simpleExit = 0.0

intold = 0.0

for x in range (minnumberofsteps, maxnumberofsteps+1):

   xtt.Start()
   simple.SetBinContent(x, simpleintegral(a, b, x, func))

   y = intold - 1.0*simpleintegral(a, b, x, func)
   y = 1.0*abs(y)/simpleintegral(a, b, x, func)
   
   simpleerr.SetBinContent(x, y)
   intold = simpleintegral(a, b, x, func)
   xtt.Stop()

   timerab = timerab + xtt.CpuTime()
   timerabr = timerabr + xtt.RealTime()
   timerSimple.SetBinContent(x, timerab) 

   flag = 0
   if y< trsh:
        flag = 1
        
   if flag == 1 and simpleExit ==0.0 and x>2:
        simpleExit = x
        print "Simple done: ",simpleExit
        print " At value: ",y
     

xac = TStopwatch()
timerac = 0.0
timeracr = 0.0
DoneTrap = 0
trapold = 0.0
trapezoidtimer = TH1F("traptimer", "traptimer", maxnumberofsteps, minnumberofsteps,maxnumberofsteps)
trapezoid = TH1F("trap", "trap", maxnumberofsteps, minnumberofsteps,maxnumberofsteps)
trapezoiderr = TH1F("traperr", "traperr", maxnumberofsteps, minnumberofsteps,maxnumberofsteps)
for x in range (minnumberofsteps, maxnumberofsteps+1):

    xac.Start()
    trapezoid.SetBinContent(x, trapezoidalint(a, b, x, func))
    y = 1.0*trapezoidalint(a, b, x, func) - trapold
    y = 1.0*abs(y)/trapezoidalint(a, b, x, func)
    trapezoiderr.SetBinContent(x,y)
    xac.Stop()

    trapold = trapezoidalint(a, b, x, func)
                               
    timerac = timerac + xac.CpuTime()
    timeracr = timeracr + xac.RealTime()
    trapezoidtimer.SetBinContent(x, timerac)

    flag = 0
    if abs(y)< trsh:
        flag = 1


    if flag == 1 and DoneTrap == 0.0 and x>1:
        DoneTrap = x
        print "DoneTrap: ",DoneTrap
        print " Error Value of Trap: ",y
        
         



timerSimple.SetLineColor(kBlue)
timerSimple.SetLineWidth(2)
trapezoidtimer.SetLineColor(kRed)
trapezoidtimer.SetLineWidth(2)

simple.SetLineColor(kBlue)
simple.SetLineWidth(2)
simpleerr.SetLineColor(kBlue)
simpleerr.SetLineWidth(2)

trapezoid.SetLineColor(kRed)
trapezoid.SetLineWidth(2)
trapezoiderr.SetLineColor(kRed)
trapezoiderr.SetLineWidth(2)





binbreak = 0





# mcav will store the average of the function after n steps, this is Fk
# mcrms will store the RMS of (Fk)
# mc is just the value of monte carlo without averaging. 


mcRMSSafonov = TH1F("mcRMSSafonov", "mcRMSSafonov", maxnumberofsteps, minnumberofsteps,maxnumberofsteps)

mc = TH1F("mc", "mc", maxnumberofsteps, minnumberofsteps,maxnumberofsteps)
mcav = TH1F("monte", "monte", maxnumberofsteps, minnumberofsteps,maxnumberofsteps)
mcRmsErr = TH1F("monteerr", "monteerr", maxnumberofsteps, minnumberofsteps,maxnumberofsteps)
mcrms = TH1F("mcrms", "mcrms", maxnumberofsteps, minnumberofsteps,maxnumberofsteps)


mctimer = TH1F("mctimer", "mctimer", maxnumberofsteps, minnumberofsteps, maxnumberofsteps)

xmct = TStopwatch()
timermc = 0.0
timermcr = 0.0


mc.SetLineColor(kOrange+4)
mc.SetLineWidth(2)
mcav.SetLineColor(kBlack)
mcrms.SetLineColor(kBlack)
mcrms.SetLineWidth(2)
mcav.SetLineWidth(2)
mcRMSSafonov.SetLineColor(kCyan+3)
mcRMSSafonov.SetLineWidth(2)
mcRmsErr.SetLineColor(kGreen+2)
mcRmsErr.SetLineWidth(2)


MCDonee = 0
temp = 0.0
temperr = 0.0
avMc = 0.0
# This works for normal MC withouto weight
for x1 in range (minnumberofsteps, maxnumberofsteps+1):

    xmct.Start()
    evalx = 1.0*MCSimpleIntegral(a, b, x1, func, 0)
    mc.SetBinContent(x1, evalx)

    
    temp = (temp + 1.0*evalx*evalx)
    mcrms.SetBinContent(x1, sqrt(temp/x1))
    avMc = (avMc + 1.0*evalx )
    mcav.SetBinContent(x1, avMc/x1)


    
    evalxerr = 1.0*MCSimpleIntegral(a, b, x1, func, 1)          # This is Delta F ^2

    variable = (sqrt(evalxerr)/(avMc/x1))/(x1)
    mcRmsErr.SetBinContent(x1, variable)            # This stores Delta I / I


    mcRMSSafonov.SetBinContent(x1, sqrt(evalxerr))

    xmct.Stop()
    timermc = timermc + xmct.CpuTime()
    timermcr = timermcr + xmct.RealTime()

    mctimer.SetBinContent(x1, timermc)



    if variable < trsh and MCDonee == 0 and x1 >5:     
        MCDonee = x1

print "MC Done at ",MCDonee
        

mctimer.SetLineColor(kBlack)
mctimer.SetLineWidth(2)


# Do the same for weighted function.
mcw = TH1F("mcw", "mcw", maxnumberofsteps, minnumberofsteps,maxnumberofsteps)
mcavw = TH1F("montew", "montew", maxnumberofsteps, minnumberofsteps,maxnumberofsteps)
mcRmsErrw = TH1F("monteerrw", "monteerrw", maxnumberofsteps, minnumberofsteps,maxnumberofsteps)
mcrmsw = TH1F("mcrmsw", "mcrmsw", maxnumberofsteps, minnumberofsteps,maxnumberofsteps)


timermcwhh = TH1F("timermcwh", "timermcwh", maxnumberofsteps, minnumberofsteps, maxnumberofsteps)


xmctwh = TStopwatch()
timermcwh = 0.0
timermcwhr = 0.0



mcrmswSaf = TH1F("mcrmswsaf", "mcrmswsaf", maxnumberofsteps, minnumberofsteps,maxnumberofsteps)

mcavw.SetLineColor(kBlack)
mcavw.SetLineWidth(2)
mcrmsw.SetLineColor(kBlack)
mcrmswSaf.SetLineColor(kCyan+3)
mcrmswSaf.SetLineWidth(2)
mcrmsw.SetLineWidth(2)
mcRmsErrw.SetLineColor(kPink+9)
mcRmsErrw.SetLineWidth(2)
mcw.SetLineColor(kOrange+4)
mcw.SetLineWidth(2)


timermcwhh.SetLineColor(kPink+9)
timermcwhh.SetLineWidth(2)

tempw = 0.0
McAvw = 0.0
errw = 0.0

MCWDone = 0

for x2 in range (minnumberofsteps, maxnumberofsteps+1):

    xmctwh.Start()
    
    evalx2 = MCIntegralWeighted(a, b, x2, func, weight, 0)    
    mcw.SetBinContent(x2, evalx2)
    tempw = tempw + evalx2*evalx2
    mcrmsw.SetBinContent(x2, sqrt(tempw/x2))
    McAvw = McAvw + evalx2
    mcavw.SetBinContent(x2, McAvw/x2)

    
    evalerrx2 = MCIntegralWeighted(a, b, x2, func, weight, 1)


    variablex2 = (sqrt(evalx2)/(McAvw/(x2)))/(x2)
    
    mcRmsErrw.SetBinContent(x2, variablex2)
    

    mcrmswSaf.SetBinContent(x2, sqrt(evalerrx2))

    xmctwh.Stop()

    timermcwh = timermcwh + xmctwh.CpuTime()
    timermcwhr = timermcwhr + xmctwh.RealTime()
    
    timermcwhh.SetBinContent(x2, timermcwh)


    flag = 0
    if variablex2 < trsh:
        flag = 1
        
    if flag == 1 and MCWDone == 0 and x2 > 5:
        MCWDone = x2


print "MCW Done at: ",MCWDone
if binbreak ==0:
    binbreak = maxnumberofsteps




if simpleExit ==0:
    simpleExit = maxnumberofsteps
if DoneTrap == 0:
    Donetrap = maxnumberofsteps

if MCDonee == 0:
    MCDonee = maxnumberofsteps

if MCWDone == 0:
    MCWDone = maxnumberofsteps
    

# Reference value
 
expectedvalue= TH1F("expected", "expected", maxnumberofsteps, minnumberofsteps,maxnumberofsteps)
expectedvalue.SetLineColor(kGreen+2)
expectedvalue.SetLineWidth(4)

for x in range (minnumberofsteps, maxnumberofsteps+1):
    expectedvalue.SetBinContent(x, func.Integral(a, b))



miny = 0
maxy = 0
if func.Integral(a,b) > 0:
    maxy = 2.60*func.Integral(a, b)
else:
    miny = 2.60*func.Integral(a, b)


    

b1 = TH2F("b1","b1",binbreak,0,binbreak,50,miny, maxy)
text = TLatex(binbreak*0.2, maxy*0.75, " f(x) = "+ff)
text2 = TLatex(binbreak*0.2, maxy*0.62, " %s < x < %s"%(a,b))



########################## Integral Convergence Comparison ########################333


legend = TLegend(0.52,0.11,0.89,0.35)
legend.SetFillColor(ROOT.kWhite)
legend.SetMargin(0.25)
legend.SetHeader(" Integral Method:")
legend.AddEntry(simple,"Midpoint", "l")
legend.AddEntry(trapezoid,"Trapezoid","l")
legend.AddEntry(mcav,"Monte Carlo (Av.)","l")


b1.GetXaxis().SetTitle("Number of Steps")
b1.GetYaxis().SetTitle("Integral Value")
b1.SetTitle("Numerical Integration: Convergence Comparison")
b1.SetStats(0)




b1.Draw()

legend.Draw("same")

trapezoid.Draw("same")

text.Draw("same")
text2.Draw("same")
simple.Draw("same")
mcav.Draw("same")
c1.SaveAs("Lab1_IntegralComparison.pdf")




############################## Times ##################################################


nnn = timermcwhh.GetMaximum()

if mctimer.GetMaximum() > nnn:
    nnn = mctimer.GetMaximum()



textab = TLatex(binbreak*0.06, 0.62*nnn, " f(x) = "+ff)
textac = TLatex(binbreak*0.05, 0.54*nnn, " %s < x < %s"%(a,b))

b3 = TH2F("b3","b3",100,0,maxnumberofsteps,200,0, 1.1*nnn)
b3.GetXaxis().SetTitle("Number of Steps")
b3.GetYaxis().SetTitle("Cumulative CPU Time (s)")
b3.SetTitle("Numerical Integration: Time Comparison")
b3.SetStats(0)


legendt = TLegend(0.11,0.61,0.42,0.88)
legendt.SetFillColor(ROOT.kWhite)
legendt.SetMargin(0.25)
legendt.SetHeader(" Integral Method:")
legendt.AddEntry(timerSimple,"Midpoint (%d)"%simpleExit, "l")
legendt.AddEntry(trapezoidtimer,"Trapezoid (%d)"%DoneTrap,"l")
legendt.AddEntry(mctimer,"Monte Carlo Av. (%d)"%MCDonee,"l")
legendt.AddEntry(timermcwhh,"Monte Carlo W. (%d)"%MCWDone,"l")






b3.Draw()
timerSimple.Draw("same")
trapezoidtimer.Draw("same")
mctimer.Draw("same")
timermcwhh.Draw("same")
legendt.Draw("same")
textab.Draw("same")
textac.Draw("same")

c1.SaveAs("IntegrationTimer.pdf")

########################### MC Data Information ##########################################


b1.SetTitle("Monte Carlo Comparison")
b1.GetYaxis().SetTitle("Value")
b1.SetStats(0)
legend2 = TLegend(0.52,0.11,0.89,0.35)
legend2.SetFillColor(ROOT.kWhite)
legend2.SetMargin(0.25)
legend2.SetHeader("Monte Carlo Value:")
#legend2.AddEntry(mc,"Monte Carlo","l")
legend2.AddEntry(mcav,"Average","l")
legend2.AddEntry(mcRMSSafonov,"RMS(\sigma)","l")
legend2.AddEntry(mcRmsErr, "\Delta I /I ", "l")
#legend2.AddEntry(expectedvalue, "Reference Value", "l")



b1.Draw()

mcRMSSafonov.Draw("same")
mcRmsErr.Draw("same")
mcav.Draw("same")
legend2.Draw("same")
text.Draw("same")
text2.Draw("same")

c1.SaveAs("Lab1_MonteCarlo.pdf")



########################### MC Weighted  Data Information ##########################################


b1.SetTitle("Monte Carlo Weighted Comparison")
b1.GetYaxis().SetTitle("Value")
b1.SetStats(0)
legend3 = TLegend(0.52,0.11,0.89,0.35)
legend3.SetFillColor(ROOT.kWhite)
legend3.SetMargin(0.25)
legend3.SetHeader(" Monte Carlo Value: ")
#legend3.AddEntry(mcw,"Monte Carlo","l")
legend3.AddEntry(mcavw,"Average of I","l")
legend3.AddEntry(mcrmswSaf,"RMS (\sigma)","l")
legend3.AddEntry(mcRmsErrw, "\Delta I /I ", "l")
#legend3.AddEntry(expectedvalue, "Reference Value", "l")




b1.Draw()
text.Draw("same")
text2.Draw("same")
mcrmswSaf.Draw("same")
mcRmsErrw.Draw("same")
#mcw.Draw("same")
mcavw.Draw("same")
legend3.Draw("same")

c1.SaveAs("Lab1_MonteCarlo_Weighted.pdf")

############################ Error For integrals ########################33


trapezoiderr.SetTitle("Numerical Integration: \Delta I / I Comparison ")
trapezoiderr.GetXaxis().SetTitle("Number of Steps")
trapezoiderr.GetYaxis().SetTitle(" \Delta I / I ")
trapezoiderr.SetStats(0)


mcRmsErrw.SetTitle("Numerical Integration: \Delta I / I Comparison ")
mcRmsErrw.GetXaxis().SetTitle("Number of Steps")
mcRmsErrw.GetYaxis().SetTitle(" \Delta I / I ")
mcRmsErrw.SetStats(0)



c1.SetLogy()


nnn = trapezoiderr.GetMaximum()

if mcRmsErrw.GetMaximum() > nnn:
    nnn = mcRmsErrw.GetMaximum()

    mcRmsErrw.Draw()
    trapezoiderr.Draw("same")
    simpleerr.Draw("same")
    mcRmsErr.Draw("same")
    mcRmsErrw.Draw("same")



else:

    
    trapezoiderr.Draw()
    
    simpleerr.Draw("same")
    mcRmsErr.Draw("same")
    mcRmsErrw.Draw("same")
    
textab2 = TLatex(binbreak*0.6, 0.82*nnn, " f(x) = "+ff)
textac2 = TLatex(binbreak*0.6, 0.35*nnn, " %s < x < %s"%(a,b))




lege = TLegend(0.16,0.65,0.46,0.88)
lege.SetFillColor(ROOT.kWhite)
lege.SetMargin(0.25)
lege.SetHeader(" Threshold:  \Delta I/I < %.3f "%trsh)
lege.AddEntry(simpleerr, "Midpoint: %d steps"%simpleExit, "l")
lege.AddEntry(trapezoiderr, "Trapezoid: %d steps"%DoneTrap, "l")
lege.AddEntry(mcRmsErr, "Monte Carlo: %d steps"%MCDonee, "l")
lege.AddEntry(mcRmsErrw, "MC Weight.: %d steps"%MCWDone, "l")

lege.Draw("same")
textab2.Draw("same")
textac2.Draw("same")
c1.SaveAs("Error_Comparison.pdf")
