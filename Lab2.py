import sys
from ROOT import *


c1 = TCanvas()
c1.SetGridx()
c1.SetGridy()
c1.SetTickx()
c1.SetTicky()


f = TFile("Lab_Assignment_004.root")
h = f.Get("h1")

p0 = 2.5
p1 = 25.0
p2 = 20.4
p3 = 4.8

mark = TMarker(1, 40, 34)
errm = TMarker(1, 40, 28)
mk = TMarker(1, 20, 7)


thresholdchisquared = 0.5
thresholdlikelihood = 0.5
#h.Draw()

#____________________________________________-
def dop3(c0 = 0, ca= 10, m2= 1.0, doit=0):

    ac =TH1F( "ac", "ac", (ca-c0), c0/m2, ca/m2)  

    temp2 = 9999.
    binsaved2 = 0

    for x in range (c0, ca):
        a2 = 1.0*x/m2
        ff = "(%f) + %f*exp( (-1.0)*(x-%f)*(x-%f)/(2*(%f)*(%f)))"%(p0,  p1, p2, p2, a2, a2)

        func = TF1("func", ff, 0, 50)
        sum2 = 0.0
                
        for i in range (1, 51):
            sum2 = sum2 + (  (  func.Eval(i, 0, 0) - h.GetBinContent(i)  )* (  func.Eval(i, 0, 0) - h.GetBinContent(i)  ) /(2*(func.Eval(i,0,0) ) ) )

        if sum2 < abs(temp2):

            temp2 = sum2
            binsaved2 = x/m2

        ac.SetBinContent(x-c0, sum2)
        #ac.Fill(a2, sum2)

    #print "Chi square min for p3: ",temp2
    #rint " value of p3: ",binsaved2



    binerr = 0.0
    binerrup = 0.0
    flag = 0
    flags = 0
    sumaa = 0.0
    sumaaup = 0.0
    pre = "Longg"



    if doit==1:
        pre = ""
        for x in range (binsaved2*m2, ca):
            a2 = 1.0*x/m2
            ff = "(%f) + %f*exp( (-1.0)*(x-%f)*(x-%f)/(2*(%f)*(%f)))"%(p0,  p1, p2, p2, a2, a2)
            func = TF1("func", ff, 0, 50)
            sum2 = 0.0
            temporal1 = 0.0
            for i in range (1, 51):
                temporal1 = temporal1 +  ((func.Eval(i,0,0)-h.GetBinContent(i))*(func.Eval(i,0,0)-h.GetBinContent(i))/(2*(func.Eval(i,0,0))))
            if temporal1 > ( temp2 + thresholdchisquared ) and flag == 0:
                flag = 1
                print (" Break at :",a2)
                print (" Value of Chi plus: ",temporal1)
                binerrup = a2
                sumaaup = temporal1
        for x in range (c0, binsaved2*m2):
            a2 = (1.0*c0)/m2 + binsaved2 - 1.0*(x)/m2
            ff = "(%f) + %f*exp( (-1.0)*(x-%f)*(x-%f)/(2*(%f)*(%f)))"%(p0,  p1, p2, p2, a2, a2)
            func = TF1("func", ff, 0, 50)
            sum2 = 0.0
            temporal2 = 0.0     
            for i in range (1, 51):
                temporal2 = temporal2+((func.Eval(i,0,0)-h.GetBinContent(i))*(func.Eval(i,0,0)-h.GetBinContent(i))/(2*(func.Eval(i,0,0))))
            if temporal2 > (temp2 + thresholdchisquared) and flags == 0:
                flags = 1
                print (" Break at :",a2)
                print (" Value of Chi minus: ",temporal2)
                binerr = a2
                sumaa = temporal2
    ac.SetLineColor(kOrange)
    ac.SetLineWidth(3)
    ac.SetTitle("Calculating p_{3} from \Chi^{2}")
    ac.GetYaxis().SetTitle("\Chi^{2}")
    ac.GetXaxis().SetTitle(" p_{3} ")
    ac.SetStats(0)
    ac.Draw("")
    
    text = TLatex(binsaved2*0.98, 0.9*temp2, " p_{3} = %.2f "%binsaved2)
    mark.DrawMarker(binsaved2, temp2)
    text.Draw("same")
    if doit == 1:
        errm.DrawMarker(binerr, sumaa)
        errm.DrawMarker(binerrup, sumaaup)
        text2 = TLatex(binsaved2*0.65, 1.022*temp2, " \Delta p_{<} = %.2f "%abs(binerr-binsaved2))
        text3 = TLatex(binsaved2*0.65, 1.018*temp2, " \Delta p_{>} = %.2f "%abs(binerrup-binsaved2))
        text2.Draw("same")
        text3.Draw("same")
    c1.SaveAs(pre+"P3.pdf")


    temp25 = 999999.
    xinit = (ca-c0)/2       # Dummy initial point
    xtemp = 0
    if doit ==0:
        for k in range (0, 5000):
            for x in range (xinit - 10, xinit + 10):
                yy = ac.GetBinContent(x)
                if abs(yy) < abs(temp25):
                    temp25 = yy
                    xtemp = x
                    valuexx = (1.0*c0 + 1.0*x)/m2
                    print (" Local min: ",temp25)
                    print (" At position: %f"%valuexx)
            if xtemp == xinit:
                valuexx = (1.0*c0 + 1.0*xinit)/m2
                print ("Abs Min Found: ",temp25)
                print ("At position: ",valuexx)
                break
            xinit = xtemp
            
        ac.SetTitle(" Minimizing p_{3}")
        ac.Draw("")
        text.Draw("same")
        binsaved2 = 1.0*binsaved2
        mark.DrawMarker(binsaved2, temp25)
        c1.SaveAs(pre+"Min_P3.pdf")
    
#_________________________________________________________________________
def dop2(c0 = 0, ca= 10, m2= 1.0, doit=0):
    #ac =TH2F( "ac", "ac", (ca-c0), c0/m2, ca/m2, 500, 0, 50)


    ac =TH1F( "ac", "ac", (ca-c0), c0/m2, ca/m2)
    temp2 = 9999.
    binsaved2 = 0

    for x in range (c0, ca):
        a2 = 1.0*x/m2
        ff = "(%f) + %f*exp( - (x-%f)*(x-%f)/(2*(%f)*(%f)))"%(p0, p1, a2, a2, p3, p3) 

        func = TF1("func", ff, 0, 50)
        sum2 = 0.0
                
        for i in range (1, 51):
            sum2 = sum2 + (  (  func.Eval(i, 0, 0) - h.GetBinContent(i)  )* (  func.Eval(i, 0, 0) - h.GetBinContent(i)  ) /(2*(func.Eval(i,0,0) ) ) )

        if sum2 < abs(temp2):

            temp2 = sum2
            binsaved2 = x/m2

        ac.SetBinContent(x-c0, sum2)
        #ac.Fill(a2, sum2)

    print ("Chi square min for p2: ",temp2)
    print (" value of p2: ",binsaved2)


    binerr = 0.0
    binerrup = 0.0
    flag = 0
    flags = 0
    sumaa = 0.0
    sumaaup = 0.0
    pre = "Longg"

    if doit==1:
        pre = ""
        
        for x in range (binsaved2*m2, ca):
            a2 = 1.0*x/m2
            ff = "(%f) + %f*exp( - (x-%f)*(x-%f)/(2*(%f)*(%f)))"%(p0, p1, a2, a2, p3, p3)
            func = TF1("func", ff, 0, 50)
            sum2 = 0.0
            temporal1 = 0.0
            for i in range (1, 51):
                temporal1 = temporal1 +  (  (  func.Eval(i, 0, 0) - h.GetBinContent(i)  )* (  func.Eval(i, 0, 0) - h.GetBinContent(i)  ) /(2*(func.Eval(i,0,0) ) ) )
            if temporal1 > ( temp2 +thresholdchisquared ) and flag == 0:
                flag = 1
                print (" Break at :",a2)
                print (" Value of Chi plus: ",temporal1)
                binerrup = a2
                sumaaup = temporal1

        for x in range (c0, binsaved2*m2):
            a2 = (1.0*c0)/m2 + binsaved2 - 1.0*(x)/m2
            ff = "(%f) + %f*exp( - (x-%f)*(x-%f)/(2*(%f)*(%f)))"%(p0, p1, a2, a2, p3, p3)
            func = TF1("func", ff, 0, 50)
            sum2 = 0.0
            temporal2 = 0.0     
            for i in range (1, 51):
                temporal2 = temporal2+ (  (  func.Eval(i, 0, 0) - h.GetBinContent(i)  )* (  func.Eval(i, 0, 0) - h.GetBinContent(i)  ) /(2*(func.Eval(i,0,0) ) ) )
            if temporal2 > (temp2 +thresholdchisquared) and flags == 0:
                flags = 1
                print (" Break at :",a2)
                print (" Value of Chi minus: ",temporal2)
                binerr = a2
                sumaa = temporal2
      

    ac.SetLineColor(kGreen+2)
    ac.SetLineWidth(3)
    ac.SetTitle("Calculating p_{2} from \Chi^{2}")
    ac.GetYaxis().SetTitle("\Chi^{2}")
    ac.GetXaxis().SetTitle(" p_{2} ")
    ac.SetStats(0)
    ac.Draw("")
    
    text = TLatex(binsaved2*0.65, 10*temp2, " p_{2} = %.2f "%binsaved2)
    mark.DrawMarker(binsaved2, temp2)

    errm.DrawMarker(binerr, sumaa)
    errm.DrawMarker(binerrup, sumaaup)

    text.Draw("same")

    if doit == 1:
        text2 = TLatex(binsaved2*0.65, 22*temp2, " \Delta p_{<} = %.2f "%abs(binerr-binsaved2))
        text3 = TLatex(binsaved2*0.65, 18*temp2, " \Delta p_{>} = %.2f "%abs(binerrup-binsaved2))
        text2.Draw("same")
        text3.Draw("same")


    c1.SaveAs(pre+"P2.pdf")

    if doit ==0:
        f = TFile(pre+"P2.root", "new")
        ac.Write()
        f.Close()



#_________________________________________________________________________
def dop1(c0 = 0, ca= 10, m2= 1.0, doit=0):

    ac =TH1F( "ac", "ac", (ca-c0), c0/m2, ca/m2)#, 500, 0, 50)

    temp2 = 9999.
    binsaved2 = 0

    for x in range (c0, ca):
        a2 = 1.0*x/m2
        ff = "(%f) + %f*exp( - (x-%f)*(x-%f)/(2*(%f)*(%f)))"%(p0, a2, p2, p2, p3, p3) 

        func = TF1("func", ff, 0, 50)
        sum2 = 0.0
                
        for i in range (1, 51):
            sum2 = sum2 + (  (  func.Eval(i, 0, 0) - h.GetBinContent(i)  )* (  func.Eval(i, 0, 0) - h.GetBinContent(i)  ) /(2*(func.Eval(i,0,0) ) ) )

        if sum2 < abs(temp2):

            temp2 = sum2
            binsaved2 = x/m2

        ac.SetBinContent(x-c0, sum2)
        #ac.Fill(a2, sum2)

    print ("Chi square min for p1: ",temp2)
    print (" value of p1: ",binsaved2)



    binerr = 0.0
    binerrup = 0.0
    flag = 0
    flags = 0
    sumaa = 0.0
    sumaaup = 0.0
    pre = "Longg"

    if doit==1:
        pre = ""
        
        for x in range (binsaved2*m2, ca):
            a2 = 1.0*x/m2
            ff = "(%f) + %f*exp( - (x-%f)*(x-%f)/(2*(%f)*(%f)))"%(p0, a2, p2, p2, p3, p3)
            func = TF1("func", ff, 0, 50)
            sum2 = 0.0
            temporal1 = 0.0
            for i in range (1, 51):
                temporal1 = temporal1 +  (  (  func.Eval(i, 0, 0) - h.GetBinContent(i)  )* (  func.Eval(i, 0, 0) - h.GetBinContent(i)  ) /(2*(func.Eval(i,0,0) ) ) )
            if temporal1 > ( temp2 +thresholdchisquared ) and flag == 0:
                flag = 1
                print (" Break at :",a2)
                print (" Value of Chi plus: ",temporal1)
                binerrup = a2
                sumaaup = temporal1

        for x in range (c0, binsaved2*m2):
            a2 = (1.0*c0)/m2 + binsaved2 - 1.0*(x)/m2
            ff = "(%f) + %f*exp( - (x-%f)*(x-%f)/(2*(%f)*(%f)))"%(p0, a2, p2, p2, p3, p3)
            func = TF1("func", ff, 0, 50)
            sum2 = 0.0
            temporal2 = 0.0     
            for i in range (1, 51):
                temporal2 = temporal2+ (  (  func.Eval(i, 0, 0) - h.GetBinContent(i)  )* (  func.Eval(i, 0, 0) - h.GetBinContent(i)  ) /(2*(func.Eval(i,0,0) ) ) )
            if temporal2 > (temp2 +thresholdchisquared) and flags == 0:
                flags = 1
                print (" Break at :",a2)
                print (" Value of Chi minus: ",temporal2)
                binerr = a2
                sumaa = temporal2
      

    ac.SetLineColor(kBlue)
    ac.SetLineWidth(3)
    ac.SetTitle("Calculating p_{1} from \Chi^{2}")
    ac.GetYaxis().SetTitle("\Chi^{2}")
    ac.GetXaxis().SetTitle(" p_{1} ")
    ac.SetStats(0)
    ac.Draw("")
    
    text = TLatex(binsaved2*0.9, 0.4*temp2, " p_{1} = %.2f "%binsaved2)
    mark.DrawMarker(binsaved2, temp2)

    errm.DrawMarker(binerr, sumaa)
    errm.DrawMarker(binerrup, sumaaup)

    text.Draw("same")

    if doit == 1:
        text2 = TLatex(binsaved2*0.85, 2.5*temp2, " \Delta p_{<} = %.2f "%abs(binerr-binsaved2))
        text3 = TLatex(binsaved2*0.85, 1.9*temp2, " \Delta p_{>} = %.2f "%abs(binerrup-binsaved2))
        text2.Draw("same")
        text3.Draw("same")

    c1.SaveAs(pre+"P1.pdf")


    if doit ==0:
        f = TFile(pre+"P1.root", "new")
        ac.Write()
        f.Close()  

#_________________________________________________________________________
def dop0(c0 = 0, ca= 10, m2= 1.0, doit=0):

    ac =TH1F( "ac", "ac", (ca-c0), c0/m2, ca/m2)#, 500, 0, 50)

    temp2 = 9999.
    binsaved2 = 0

    for x in range (c0, ca):
        a2 = 1.0*x/m2
        ff = "(%f) + %f*exp( (-1.0)*(x-%f)*(x-%f)/(2*(%f)*(%f)))"%(a2, p1, p2, p2, p3, p3) 

        func = TF1("func", ff, 0, 50)
        sum2 = 0.0
                
        for i in range (1, 51):
            sum2 = sum2 + (  (  func.Eval(i, 0, 0) - h.GetBinContent(i)  )* (  func.Eval(i, 0, 0) - h.GetBinContent(i)  ) /(2*(func.Eval(i,0,0) ) ) )

        if sum2 < abs(temp2):

            temp2 = sum2
            binsaved2 = x/m2

        ac.SetBinContent(x-c0, sum2)
        #ac.Fill(a2, sum2)

    print ("Chi square min for p0: ",temp2)
    print (" value of p0: ",binsaved2)



    binerr = 0.0
    binerrup = 0.0
    flag = 0
    flags = 0
    sumaa = 0.0
    sumaaup = 0.0
    pre = "Longg"

    if doit==1:
        pre = ""
        
        for x in range (binsaved2*m2, ca):
            a2 = 1.0*x/m2
            ff = "(%f) + %f*exp( (-1.0)*(x-%f)*(x-%f)/(2*(%f)*(%f)))"%(a2, p1, p2, p2, p3, p3) 
            func = TF1("func", ff, 0, 50)
            sum2 = 0.0
            temporal1 = 0.0
            for i in range (1, 51):
                temporal1 = temporal1 +  (  (  func.Eval(i, 0, 0) - h.GetBinContent(i)  )* (  func.Eval(i, 0, 0) - h.GetBinContent(i)  ) /(2*(func.Eval(i,0,0) ) ) )
            if temporal1 > ( temp2 +thresholdchisquared ) and flag == 0:
                flag = 1
                print (" Break at :",a2)
                print (" Value of Chi plus: ",temporal1)
                binerrup = a2
                sumaaup = temporal1

        for x in range (c0, binsaved2*m2):
            a2 = (1.0*c0)/m2 + binsaved2 - 1.0*(x)/m2
            ff = "(%f) + %f*exp( (-1.0)*(x-%f)*(x-%f)/(2*(%f)*(%f)))"%(a2, p1, p2, p2, p3, p3) 
            func = TF1("func", ff, 0, 50)
            sum2 = 0.0
            temporal2 = 0.0     
            for i in range (1, 51):
                temporal2 = temporal2+ (  (  func.Eval(i, 0, 0) - h.GetBinContent(i)  )* (  func.Eval(i, 0, 0) - h.GetBinContent(i)  ) /(2*(func.Eval(i,0,0) ) ) )
            if temporal2 > (temp2 +thresholdchisquared) and flags == 0:
                flags = 1
                print (" Break at :",a2)
                print (" Value of Chi minus: ",temporal2)
                binerr = a2
                sumaa = temporal2
      

    ac.SetTitle("Calculating p_{0} from #Chi^{2}")
    ac.GetYaxis().SetTitle("#Chi^{2}")
    ac.GetXaxis().SetTitle(" p_{0} ")
    ac.SetStats(0)
    ac.SetLineColor(kRed)
    ac.SetLineWidth(4)
    ac.Draw("")
    
    text = TLatex(binsaved2*0.75, 4*temp2, " p_{0}^{2} = %.2f "%binsaved2)
    mark.DrawMarker(binsaved2, temp2)

    errm.DrawMarker(binerr, sumaa)
    errm.DrawMarker(binerrup, sumaaup)

    text.Draw("same")

    if doit == 1:
        text2 = TLatex(binsaved2*0.75, 15*temp2, " \Delta p_{<} = %.2f "%abs(binerr-binsaved2))
        text3 = TLatex(binsaved2*0.75, 13*temp2, " \Delta p_{>} = %.2f "%abs(binerrup-binsaved2))
        text2.Draw("same")
        text3.Draw("same")

    c1.SaveAs(pre+"P0.pdf")

    if doit ==0:
        f = TFile(pre+"P0.root", "new")
        ac.Write()
        f.Close()

# ____________________________________________________________________________
def minimizerF(func = TF1("f1", "f1", 0, 1), a = 0, b = 1, n = 1000.0):

    ab = TH1F("ab", "ab", n, a, b)
    dxy =  999.
    oldpoint = 0
    for x in range (0, n):
        x2 = 1.0*a + (1.0*b - 1.0*a)*x/(1.0*n)
        y = func.Eval(x2, 0, 0)
        ab.SetBinContent(x, y)
        deriv = 1.0*func.Eval(1.0*a + (1.0*b - 1.0*a)*(1+x)/(1.0*n), 0, 0)-y
        deriv = deriv*1.0*n/(1.0*b - 1.0*a)
        if dxy/deriv < 0:
            ismin = 0
            ismax = 0
            for xy in range (x - 10, x):
                xx3 = 1.0*a + (1.0*b - 1.0*a)*xy/(1.0*n)
                yy = func.Eval(xx3, 0, 0)
                if yy < y:
                    ismax = 1
                if yy > y:
                    ismin = 1
            for xy in range (x, x + 10):
                xx3 = 1.0*a + (1.0*b - 1.0*a)*xy/(1.0*n)
                yy = func.Eval(xx3, 0, 0)
                if yy < y:
                    ismax = ismax + 1
                if yy > y:
                    ismin = ismin + 1
            exx = " ? "
            if ismin > 1:
                exx = "Minimum"
            if ismax > 1:
                exx = "Maximum"
            print ("Extreme Point: "+exx)
            print (" At x: ",x2)
            print (" Value: ",func.Eval(x2, 0, 0))
            oldpoint = func.Eval(x2,0,0)
        dxy = func.Derivative(x2)
    print ("Doesn't print for extremes at the edges of the interval") 
    ab.Draw()
    c1.SaveAs("min.pdf")





#________________________________________________________________________

def minimize(name = "p3"):

    if name == "p3":
        ac = dop3(4000, 6000, 1000., 1)

    n = ac.GetXaxis().GetNbins()
    a = ac.GetXaxis().GetXmin()


    stopsd = 0
    stopsg = 0
    minvalue = 0
    maxvalue = 0


    ada = 0.0
    
    y = ac.GetBinContent(1)

    if ac.GetBinContent(2)> y:
        for x in range (2, n+1):
            if (ac.GetBinContent(x) > ac.GetBinContent(x+1)):
                   print (" Stops growing at ",x)
                   stopsg = x


        maxvalue = ac.GetBinContent(x)

        
        for xx in range (x , x+10):

                if x + 1> n:
                    break
                
                if ac.GetBinContent(xx)<maxvalue:
                        maxvalue = maxvalue

                else:
                        maxvalue = ac.GetBinContent(xx)


    if ac.GetBinContent(2) < y:
        for x in range (2, n+1):
            
            if ac.GetBinContent(x) < ac.GetBinContent(x+1):
                #print "Stops decreasing at ",1.0*a+ (1.0*x)/n
                stopsd = 1.0*a+ (1.0*x)/n
                ada= x
                


        
        for xx in range (x, x+10):
            if ac.GetBinContent(xx)<minvalue:
               minvalue = ac.GetBinContent(xx)
               stopsd = 1.0*a + (1.0*xx)/n

    print (stopsd)
    minvalue = ac.GetBinContent(ada)
    print ("Min Value: ",minvalue)

    ac.Draw()
    text =  TLatex(stopsd, 0.99*minvalue, "X")

    text2 = TLatex(stopsd, 0.79*minvalue, "Min p_{3}: %.2f"%stopsd)

    ac.Draw()
    text.Draw("same")
    text2.Draw("same")

    c1.SaveAs("minimize_formal_"+name+".pdf")
    
    
#____________________________________________________________________________
def minimizea(name= "p3"):

    if name == "p3":
        ac = dop3(4000, 6000, 1000., 1)

    if name == "p2":
        ac = dop2(1700, 2500, 100., 1)
    n = ac.GetXaxis().GetNbins()
    a = ac.GetXaxis().GetXmin()
    maxvalue = 0.
    minvalue = 99999999.
    maxx = 0
    minx = 0 
    for x in range (0, n):
        y = ac.GetBinContent(x)
        if  y > maxvalue:
            maxvalue = y
            maxx = x

        if (abs(y) < abs(minvalue)):
            minvalue = y
            minx = x


    minx = a+(1.0*x)/n

    print ("Miny: ",minvalue)
    print ("Minx: ",minx)


    
    ac.SetLineColor(kRed)
    ac.SetLineWidth(2)
    ac.SetTitle("Minimizing "+name)
    ac.GetYaxis().SetTitle("\Chi^{2}")
    ac.GetXaxis().SetTitle(name)
    ac.SetStats(0) 
    text =  TLatex(minx, 0.99*minvalue, "X")

    text2 = TLatex(minx, 0.79*minvalue, "Min p_{3}: %.2f"%minx)

    ac.Draw()
    text.Draw("same")
    text2.Draw("same")
    
    c1.SaveAs("min_"+name+".pdf")
    

# ____________________________________________________________________________
def minimizerH1D(func = TH1F("f1", "f1", 100, 0, 1), a = 0, b = 1, n = 1000.0, name="name"):

    ab = TH1F("ab", "ab", n, a, b)
    dxy =  999.
    oldpoint = 0

    
    for x in range (0, n):
        x2 = 1.0*a + (1.0*b - 1.0*a)*x/(1.0*n)
        y = func.GetBinContent(x)
        ab.SetBinContent(x, y)

        

        deriv = 1.0*func.GetBinContent(x+1)-y
        deriv = deriv*1.0*n/(1.0*b - 1.0*a)


        
        if dxy/deriv < 0:

            ismin = 0
            ismax = 0
        
            for xy in range (x - 10, x):
                xx3 = 1.0*a + (1.0*b - 1.0*a)*xy/(1.0*n)
                yy = func.GetBinContent(xy)

                if yy < y:
                    ismax = 1

                if yy > y:
                    ismin = 1

            for xy in range (x, x + 10):
                xx3 = 1.0*a + (1.0*b - 1.0*a)*xy/(1.0*n)
                yy = func.GetBinContent(xy)

                if yy < y:
                    ismax = ismax + 1

                if yy > y:
                    ismin = ismin + 1


            exx = " ? "
            
            if ismin > 1:
                exx = "Minimum"

            if ismax > 1:
                exx = "Maximum"


            oldpoint = func.GetBinContent(x)
            print ("Extreme Point: "+exx)
            print (" At x: ",x2)
            print (" Value: ",oldpoint)

            
            
        dxy = deriv

    print ("Doesn't print for extremes at the edges of the interval" )



    
    ab.Draw()
    c1.SaveAs("min+"+name+".pdf")



#_________________________________________________________________________
def dop2p3(c0 = 0, ca= 10, m2= 1.0, a0= 0, a1 = 10, m3 = 1.0, doit = 0):

    ak =TH2F( "ak", "ak", ca-c0, c0/m2, ca/m2, a1-a0, a0/m3, a1/m3)

    temp2 = 9999.
    binsaved2 = 0.0
    binsaved3 = 0.0

    
    for x2 in range (c0, ca):
        b0 = 1.0*x2/m2

 
        for y in range (a0, a1):
            a2 = 1.0*y/m3
            ff = "(%f) + %f*exp( (-1.0)*(x-%f)*(x-%f)/(2*(%f)*(%f)))"%(p0, p1, b0, b0, a2, a2) 

            func = TF1("func", ff, 0, 41)
            sum2 = 0.0
            
            for i in range (1, 41):
                sum2 = sum2 + (  (  func.Eval(i, 0, 0) - h.GetBinContent(i)  )* (  func.Eval(i, 0, 0) - h.GetBinContent(i)  ) /(2*(func.Eval(i,0,0) ) ) )

            if sum2 < abs(temp2):

                temp2 = sum2
                binsaved2 = b0
                binsaved3 = a2


            ak.SetBinContent(x2 - c0, y-a0, sum2)

    #print "Chi square min for p2: ",temp2
    #print " value of p2: ",binsaved2
    #rint " Value of p3: ",binsaved3

    pre = "Long"
    scale1= 1.3
    scale2 = 0.7

    ak.SetTitle("Calculating p_{3} and p_{2} from #chi^{2}")
    ak.GetYaxis().SetTitle(" p_{3}")
    ak.GetXaxis().SetTitle("  p_{2} ")
    ak.SetStats(0)

    
    ak.Draw("colz")   
    mark.DrawMarker(binsaved2, binsaved3)
 


    if doit==1:
        scale1= 0.9
        scale2=0.88
        pre = ""
        for x2 in range (c0, ca):
            b0 = 1.0*x2/m2
            for y in range (a0, a1):
                a2 = 1.0*y/m3
                ff ="(%f) + %f*exp( (-1.0)*(x-%f)*(x-%f)/(2*(%f)*(%f)))"%(p0, p1, b0, b0, a2, a2) 
                func = TF1("func", ff, 0, 41)
                sum2 = 0.0
                for i in range (1, 41):
                    sum2 = sum2 +((func.Eval(i,0,0)-h.GetBinContent(i))*(func.Eval(i,0,0)-h.GetBinContent(i))/(2*(func.Eval(i,0,0))))

                if abs(sum2) < abs(temp2)+thresholdchisquared:
                    if abs(sum2) > abs(temp2 + 0.93*thresholdchisquared):
                        mk.DrawMarker(b0, a2)


    text2 = TLatex(0.96*binsaved2, binsaved3*scale1, " P2: %.2f"%binsaved2)
    text3 = TLatex(0.96*binsaved2, binsaved3*scale2, " P3: %.2f"%binsaved3)    
    text2.Draw("same")
    text3.Draw("same") 
    c1.SaveAs(pre+"P2P3.pdf")



    if doit == 0:                              
        xinit = (ca - c0)/2                     # Dummy initial parameters
        yinit = (a1 - a0)/2     
        tempp = 10000                           # Dummy initial values
        xtemp = 0
        ytemp = 0
        
        for k in range (0, 5000):
            for x2 in range (xinit - 10, xinit + 10):
                for y2 in range (yinit - 10, yinit + 10):
                    ff = ak.GetBinContent(ak.GetBin(x2, y2))
                                            
                    if ff < tempp:
                        tempp = ff          # Temporary storing a possible min
                        xtemp = x2          # And its coordinates
                        ytemp = y2

            #print "Min Local Value:",tempp
            #print "Coordinates of it: %f, %f"%(1.0*xtemp/m2, 1.0*ytemp/m2)

            # The following conditions claims that the min hasn't changed.
            if xinit == xtemp and yinit == ytemp:
                #print "Found the Absolute min at: %f, %f"%(xinit/m2, yinit/m2)
                break
            
            xinit = xtemp
            yinit = ytemp

        ak.SetTitle("Minimization of \chi^{2} as a function of p_{3} and p_{2}")
        ak.GetYaxis().SetTitle(" p_{3}")
        ak.GetXaxis().SetTitle("  p_{2} ")
        ak.SetStats(0)
        ak.Draw("colz")
        mark.DrawMarker(xinit/m2, yinit/m2)
        valueasdf = 1.0*xinit/m2
        valueasdg = 1.0*yinit/m2
        text2 = TLatex(0.96*xinit/m2, yinit*1.3/m2, " P2: %.2f"%valueasdf)
        text3 = TLatex(0.96*xinit/m2, yinit*0.7/m2, " P3: %.2f"%valueasdg)    
        text2.Draw("same")
        text3.Draw("same") 

        c1.SaveAs(pre+"Min_P2P3.pdf")




        
#_________________________________________________________________________
def dop1p2(c0 = 0, ca= 10, m2= 1.0, a0= 0, a1 = 10, m3 = 1.0, doit = 0):



    ak =TH2F( "ak", "ak", ca-c0, c0/m2, ca/m2, a1-a0, a0/m3, a1/m3)

    temp2 = 9999.
    binsaved2 = 0.0
    binsaved3 = 0.0

    
    for x2 in range (c0, ca):
        b0 = 1.0*x2/m2

 
        for y in range (a0, a1):
            a2 = 1.0*y/m3
            ff = "(%f) + %f*exp( (-1)*(x-%f)*(x-%f)/(2*(%f)*(%f)))"%(p0, b0, a2, a2, p3, p3) 

            func = TF1("func", ff, 0, 41)
            sum2 = 0.0
            
            for i in range (1, 41):
                sum2 = sum2 + (  (  func.Eval(i, 0, 0) - h.GetBinContent(i)  )* (  func.Eval(i, 0, 0) - h.GetBinContent(i)  ) /(2*(func.Eval(i,0,0) ) ) )

            if sum2 < abs(temp2):

                temp2 = sum2
                binsaved2 = b0
                binsaved3 = a2


            ak.SetBinContent(x2 - c0, y-a0, sum2) 


    #print "Chi square min for p2: ",temp2
    #print " value of p1: ",binsaved2
    #print " Value of p2: ",binsaved3


    pre = "Long"
    scale1= 1.3
    scale2 = 0.7

 
    ak.SetTitle("Calculating p_{1} and p_{2} from #chi^{2}")
    ak.GetYaxis().SetTitle(" p_{2}")
    ak.GetXaxis().SetTitle("  p_{1} ")
    ak.SetStats(0)

    ak.Draw("colz")   

    if doit==1:
        scale1=  0.97
        scale2=0.96
        pre = ""
        for x2 in range (c0, ca):
            b0 = 1.0*x2/m2
            for y in range (a0, a1):
                a2 = 1.0*y/m3
                ff ="(%f) + %f*exp( (-1)*(x-%f)*(x-%f)/(2*(%f)*(%f)))"%(p0, b0, a2, a2, p3, p3) 
                func = TF1("func", ff, 0, 41)
                sum2 = 0.0
                for i in range (1, 41):
                    sum2 = sum2 + (  (  func.Eval(i, 0, 0) - h.GetBinContent(i)  )* (  func.Eval(i, 0, 0) - h.GetBinContent(i)  ) /(2*(func.Eval(i,0,0) ) ) )

                if abs(sum2) < abs(temp2)+thresholdchisquared:
                    if abs(sum2) > abs(temp2 + 0.93*thresholdchisquared):
                        mk.DrawMarker(b0, a2)

    text2 = TLatex(0.99*binsaved2, binsaved3*scale1, " P1: %.2f"%binsaved2)
    text3 = TLatex(0.99*binsaved2, binsaved3*scale2, " P2: %.2f"%binsaved3)
    
    mark.DrawMarker(binsaved2, binsaved3)
    text2.Draw("same")
    text3.Draw("same")
    
    c1.SaveAs(pre+"P1P2.pdf")


    
    if doit == 0:                              
        xinit = (ca - c0)/2                     # Dummy initial parameters
        yinit = (a1 - a0)/2     
        tempp = 10000                           # Dummy initial values
        xtemp = 0
        ytemp = 0
        
        for k in range (0, 5000):
            for x2 in range (xinit - 10, xinit + 10):
                for y2 in range (yinit - 10, yinit + 10):
                    ff = ak.GetBinContent(ak.GetBin(x2, y2))
                                            
                    if ff < tempp:
                        tempp = ff          # Temporary storing a possible min
                        xtemp = x2          # And its coordinates
                        ytemp = y2

            #print "Min Local Value:",tempp
            #print "Coordinates of it: %f, %f"%(1.0*xtemp/m2, 1.0*ytemp/m2)

            # The following conditions claims that the min hasn't changed.
            if xinit == xtemp and yinit == ytemp:
                #print "Found the Absolute min at: %f, %f"%(xinit/m2, yinit/m2)
                break
            
            xinit = xtemp
            yinit = ytemp

        ak.SetTitle("Minimization of \chi^{2} as a function of p_{1} and p_{2}")
        ak.GetYaxis().SetTitle(" p_{2}")
        ak.GetXaxis().SetTitle("  p_{1} ")
        ak.SetStats(0)
        ak.Draw("colz")
        mark.DrawMarker(xinit/m2, yinit/m2)

        binsaved2 = 1.0*xinit/m2
        binsaved3 = 1.0*xinit/m2
        
        text2 = TLatex(0.96*xinit/m2, yinit*1.3/m2, " P1: %.2f"%binsaved2)
        text3 = TLatex(0.96*xinit/m2, yinit*0.7/m2, " P2: %.2f"%binsaved3)    
        text2.Draw("same")
        text3.Draw("same") 

        c1.SaveAs(pre+"Min_P1P2.pdf")

        

#_________________________________________________________________________________________________________________

def dop1p3(c0 = 0, ca= 10, m2= 1.0, a0= 0, a1 = 10, m3 = 1.0, doit=0):



    ak =TH2F( "ak", "ak", ca-c0, c0/m2, ca/m2, a1-a0, a0/m3, a1/m3)

    temp2 = 9999.
    binsaved2 = 0.0
    binsaved3 = 0.0

    
    for x2 in range (c0, ca):
        b0 = 1.0*x2/m2

 
        for y in range (a0, a1):
            a2 = 1.0*y/m3
            ff = "(%f) + %f*exp( (-1)*(x-%f)*(x-%f)/(2*(%f)*(%f)))"%(p0, b0, p2, p2, a2, a2) 

            func = TF1("func", ff, 0, 41)
            sum2 = 0.0
            
            for i in range (1, 41):
                sum2 = sum2 + (  (  func.Eval(i, 0, 0) - h.GetBinContent(i)  )* (  func.Eval(i, 0, 0) - h.GetBinContent(i)  ) /(2*(func.Eval(i,0,0) ) ) )

            if sum2 < abs(temp2):

                temp2 = sum2
                binsaved2 = b0
                binsaved3 = a2


            ak.SetBinContent(x2 - c0, y-a0, sum2) 


    #print "Chi square min for p2: ",temp2
    #print " value of p1: ",binsaved2
    #print " Value of p3: ",binsaved3

    pre = "Long"
    scale1= 1.3
    scale2 = 0.7
    ak.SetTitle("Calculating p_{1} and p_{3} from #chi^{2}")
    ak.GetYaxis().SetTitle(" p_{3}")
    ak.GetXaxis().SetTitle("  p_{1} ")
    ak.SetStats(0)

    ak.Draw("colz")   


    if doit==1:
        scale1= 0.9
        scale2=0.88
        pre = ""
        for x2 in range (c0, ca):
            b0 = 1.0*x2/m2
            for y in range (a0, a1):
                a2 = 1.0*y/m3
                ff ="(%f) + %f*exp( (-1)*(x-%f)*(x-%f)/(2*(%f)*(%f)))"%(p0, b0, p2, p2, a2, a2) 
                func = TF1("func", ff, 0, 41)
                sum2 = 0.0
                for i in range (1, 41):
                    sum2 = sum2 + (  (  func.Eval(i, 0, 0) - h.GetBinContent(i)  )* (  func.Eval(i, 0, 0) - h.GetBinContent(i)  ) /(2*(func.Eval(i,0,0) ) ) )

                if abs(sum2) < abs(temp2)+thresholdchisquared:
                    if abs(sum2) > abs(temp2 + 0.93*thresholdchisquared):
                        mk.DrawMarker(b0, a2)

    text2 = TLatex(0.96*binsaved2, binsaved3*scale1, " P1: %.2f"%binsaved2)
    text3 = TLatex(0.96*binsaved2, binsaved3*scale2, " P2: %.2f"%binsaved3)
    
    mark.DrawMarker(binsaved2, binsaved3)
    text2.Draw("same")
    text3.Draw("same")
    


    c1.SaveAs(pre+"P1P3.pdf")


#_________________________________________________________________________________________________________________

def dop0p1(c0 = 0, ca= 10, m2= 1.0, a0= 0, a1 = 10, m3 = 1.0, doit=0):

    ak =TH2F( "ak", "ak", ca-c0, c0/m2, ca/m2, a1-a0, a0/m3, a1/m3)

    temp2 = 9999.
    binsaved2 = 0.0
    binsaved3 = 0.0

    
    for x2 in range (c0, ca):
        b0 = 1.0*x2/m2

 
        for y in range (a0, a1):
            a2 = 1.0*y/m3
            ff ="(%f) + %f*exp( (-1)*(x-%f)*(x-%f)/(2*(%f)*(%f)))"%(b0, a2, p2, p2, p3, p3)  

            func = TF1("func", ff, 0, 41)
            sum2 = 0.0
            
            for i in range (1, 41):
                sum2 = sum2 + (  (  func.Eval(i, 0, 0) - h.GetBinContent(i)  )* (  func.Eval(i, 0, 0) - h.GetBinContent(i)  ) /(2*(func.Eval(i,0,0) ) ) )

            if sum2 < abs(temp2):

                temp2 = sum2
                binsaved2 = b0
                binsaved3 = a2


            ak.SetBinContent(x2 - c0, y-a0, sum2) 


    #print "Chi square min for p2: ",temp2
    #print " value of p0: ",binsaved2
    #print " Value of p1: ",binsaved3

    pre = "Long"
    scale1= 1.3
    scale2 = 0.7
    ak.SetTitle("Calculating p_{0} and p_{1} from #chi^{2}")
    ak.GetYaxis().SetTitle(" p_{1}")
    ak.GetXaxis().SetTitle("  p_{0}^{2} ")
    ak.SetStats(0)

    ak.Draw("colz")   


    if doit==1:
        scale1= 1.
        scale2=0.96
        pre = ""
        for x2 in range (c0, ca):
            b0 = 1.0*x2/m2
            for y in range (a0, a1):
                a2 = 1.0*y/m3
                ff ="(%f) + %f*exp( (-1)*(x-%f)*(x-%f)/(2*(%f)*(%f)))"%(b0, a2, p2, p2, p3, p3)  
                func = TF1("func", ff, 0, 41)
                sum2 = 0.0
                for i in range (1, 41):
                    sum2 = sum2 + (  (  func.Eval(i, 0, 0) - h.GetBinContent(i)  )* (  func.Eval(i, 0, 0) - h.GetBinContent(i)  ) /(2*(func.Eval(i,0,0) ) ) )

                if abs(sum2) < abs(temp2)+thresholdchisquared:
                    if abs(sum2) > abs(temp2 + 0.93*thresholdchisquared):
                        mk.DrawMarker(b0, a2)

    text2 = TLatex(1.2*binsaved2, binsaved3*scale1, " P0^{2}: %.2f"%binsaved2)
    text3 = TLatex(1.2*binsaved2, binsaved3*scale2, " P1: %.2f"%binsaved3)
    
    mark.DrawMarker(binsaved2, binsaved3)
    text2.Draw("same")
    text3.Draw("same")
    


    c1.SaveAs(pre+"P0P1.pdf")







 


#___________________________________________________
def Factorial( n = 1 ):
    if (n < 0): return 1
    x= 1.0
    for b in range (1, n+1):
        x = x*b
    return x
        
#____________________________________________-
def lp3(c0 = 0, ca= 10, m2= 1.0, doit=0):

    ac =TH1F( "ac", "ac", (ca-c0), c0/m2, ca/m2)#, 500, -200, 0)
    temp2 = -9999.
    binsaved2 = 0
    for x in range (c0, ca):
        a2 = 1.0*x/m2
        ff = "(%f) + %f*exp( (-1.0)*(x-%f)*(x-%f)/(2*(%f)*(%f)))"%(p0,  p1, p2, p2, a2, a2)
        func = TF1("func", ff, 0, 50)
        sum2 = 0.0      
        for i in range (1, 51):
            sum2 = sum2 + (h.GetBinContent(i))*log(func.Eval(i, 0, 0)) - func.Eval(i, 0, 0) - log( Factorial( h.GetBinContent(i)))
        if abs(sum2) < abs(temp2):
            temp2 = sum2
            binsaved2 = x/m2
        ac.SetBinContent(x-c0, sum2)
        #ac.Fill(a2, sum2)
    #print "Chi square min for p3: ",temp2
    #print " value of p3: ",binsaved2
    
    binerr = 0.0
    binerrup = 0.0
    flag = 0
    flags = 0

    sumaa = 0.0
    sumaup = 0.0
    pre = "L"


    if doit==1:
        pre = ""
        
        for x in range (binsaved2*m2, ca):
            a2 = 1.0*x/m2
            ff = "(%f) + %f*exp( (-1.0)*(x-%f)*(x-%f)/(2*(%f)*(%f)))"%(p0,  p1, p2, p2, a2, a2)
            func = TF1("func", ff, 0, 50)
            sum2 = 0.0
            temporal1 = 0.0
            for i in range (1, 51):
                temporal1 = temporal1 + (h.GetBinContent(i))*log(func.Eval(i, 0, 0)) - func.Eval(i, 0, 0)*log( 2.718 ) - log( Factorial( h.GetBinContent(i)))
            if temporal1 < ( temp2 -thresholdlikelihood ) and flag == 0:
                flag = 1
                #print " Break at :",a2
                #print " Value of Chi plus: ",temporal1
                binerrup = a2
                sumaaup = temporal1

        for x in range (c0, binsaved2*m2):
            a2 = (1.0*c0)/m2 + binsaved2 - 1.0*(x)/m2
            ff = "(%f) + %f*exp( (-1.0)*(x-%f)*(x-%f)/(2*(%f)*(%f)))"%(p0,  p1, p2, p2, a2, a2) 
            func = TF1("func", ff, 0, 50)
            sum2 = 0.0
            temporal2 = 0.0     
            for i in range (1, 51):
                temporal2 = temporal2 + (h.GetBinContent(i))*log(func.Eval(i, 0, 0)) - func.Eval(i, 0, 0)*log( 2.718 ) - log( Factorial( h.GetBinContent(i)))
            if temporal2 < (temp2 -thresholdlikelihood) and flags == 0:
                flags = 1
                #print " Break at :",a2
                #print " Value of Chi minus: ",temporal2
                binerr = a2
                sumaa = temporal2
      

    ac.SetTitle("Calculating p_{3} from Likelihold")
    ac.GetYaxis().SetTitle("Likelihold")
    ac.GetXaxis().SetTitle(" p_{3} ")
    ac.SetStats(0)
    ac.SetLineColor(kOrange)
    ac.SetLineWidth(4)
    ac.Draw("")
    text = TLatex(0.95*binsaved2, 0.72*temp2, " p_{3} = %.2f "%binsaved2)
    
    mark.DrawMarker(binsaved2, temp2)
    text.Draw("same")


    errm.DrawMarker(binerr, sumaa)
    errm.DrawMarker(binerrup, sumaaup)

    #print sumaa
    #print sumaaup
    if doit == 1:
        text2 = TLatex(binsaved2*4.75, 2.1*temp2, " \Delta p_{<} = %.2f "%abs(binerr-binsaved2))
        text3 = TLatex(binsaved2*4.75, 2.4*temp2, " \Delta p_{>} = %.2f "%abs(binerrup-binsaved2))
        text2.Draw("same")
        text3.Draw("same")

    c1.SaveAs(pre+"LP3.pdf")

    
#_________________________________________________________________________
def lp2(c0 = 0, ca= 10, m2= 1.0, doit=0):

    ac =TH1F( "ac", "ac", (ca-c0), c0/m2, ca/m2)#, 500, -200, 0)
    temp2 = -9999.
    binsaved2 = 0
    for x in range (c0, ca):
        a2 = 1.0*x/m2
        ff = "(%f) + %f*exp( - (x-%f)*(x-%f)/(2*(%f)*(%f)))"%(p0, p1, a2, a2, p3, p3)
        func = TF1("func", ff, 0, 50)
        sum2 = 0.0      
        for i in range (1, 51):
            sum2 = sum2 + (h.GetBinContent(i))*log(func.Eval(i, 0, 0)) - func.Eval(i, 0, 0) - log( Factorial( h.GetBinContent(i)))
        if abs(sum2) < abs(temp2):
            temp2 = sum2
            binsaved2 = x/m2
        ac.SetBinContent(x-c0, sum2)
        #ac.Fill(a2, sum2)
    #print "Chi square min for p2: ",temp2
    #print " value of p2: ",binsaved2
    
    binerr = 0.0
    binerrup = 0.0
    flag = 0
    flags = 0

    sumaa = 0.0
    sumaup = 0.0
    
    pre = "L"
    
    if doit==1:
        pre = ""
        for x in range (binsaved2*m2, ca):
            a2 = 1.0*x/m2
            ff = "(%f) + %f*exp( - (x-%f)*(x-%f)/(2*(%f)*(%f)))"%(p0, p1, a2, a2, p3, p3)
            func = TF1("func", ff, 0, 50)
            sum2 = 0.0
            temporal1 = 0.0
            for i in range (1, 51):
                temporal1 = temporal1 + (h.GetBinContent(i))*log(func.Eval(i, 0, 0)) - func.Eval(i, 0, 0)*log( 2.718 ) - log( Factorial( h.GetBinContent(i)))
            if temporal1 < ( temp2 -thresholdlikelihood ) and flag == 0:
                flag = 1
                #print " Break at :",a2
                #print " Value of Chi plus: ",temporal1
                binerrup = a2
                sumaaup = temporal1

        for x in range (c0, binsaved2*m2):
            a2 = (1.0*c0)/m2 + binsaved2 - 1.0*(x)/m2
            ff  = "(%f) + %f*exp( - (x-%f)*(x-%f)/(2*(%f)*(%f)))"%(p0, p1, a2, a2, p3, p3)
            func = TF1("func", ff, 0, 50)
            sum2 = 0.0
            temporal2 = 0.0     
            for i in range (1, 51):
                temporal2 = temporal2 + (h.GetBinContent(i))*log(func.Eval(i, 0, 0)) - func.Eval(i, 0, 0)*log( 2.718 ) - log( Factorial( h.GetBinContent(i)))
            if temporal2 < (temp2 -thresholdlikelihood) and flags == 0:
                flags = 1
                #print " Break at :",a2
                #print " Value of Chi minus: ",temporal2
                binerr = a2
                sumaa = temporal2
  

    ac.SetTitle("Calculating p_{2} from Likelihold")
    ac.GetYaxis().SetTitle("Likelihold")
    ac.GetXaxis().SetTitle(" p_{2} ")
    ac.SetStats(0)
    ac.SetLineColor(kGreen+2)
    ac.SetLineWidth(4)
    
    text = TLatex(0.75*binsaved2, 0.76*temp2, " p_{2} = %.2f "%binsaved2)
    
    ac.Draw("")
    mark.DrawMarker(binsaved2, temp2)
    text.Draw("same")


    errm.DrawMarker(binerr, sumaa)
    errm.DrawMarker(binerrup, sumaaup)
                    
    if doit == 1:
        text2 = TLatex(binsaved2*0.75, 2.9*temp2, " \Delta p_{<} = %.2f "%abs(binerr-binsaved2))
        text3 = TLatex(binsaved2*0.75, 3.3*temp2, " \Delta p_{>} = %.2f "%abs(binerrup-binsaved2))
        text2.Draw("same")
        text3.Draw("same")


        
    c1.SaveAs(pre+"LP2.pdf")




#_________________________________________________________________________
def lp1(c0 = 0, ca= 10, m2= 1.0, doit=0):

    ac =TH1F( "ac", "ac", (ca-c0), c0/m2, ca/m2)#, 500, -200, 0)
    temp2 = -9999.
    binsaved2 = 0
    for x in range (c0, ca):
        a2 = 1.0*x/m2
        ff = "(%f) + %f*exp( - (x-%f)*(x-%f)/(2*(%f)*(%f)))"%(p0, a2, p2, p2, p3, p3)
        func = TF1("func", ff, 0, 50)
        sum2 = 0.0      
        for i in range (1, 51):
            sum2 = sum2 + (h.GetBinContent(i))*log(func.Eval(i, 0, 0)) - func.Eval(i, 0, 0) - log( Factorial( h.GetBinContent(i)))
        if abs(sum2) < abs(temp2):
            temp2 = sum2
            binsaved2 = x/m2
        ac.SetBinContent(x-c0, sum2)
        #ac.Fill(a2, sum2)
    #print "Chi square min for p1: ",temp2
    #print " value of p1: ",binsaved2
    
    binerr = 0.0
    binerrup = 0.0
    flag = 0
    flags = 0

    sumaa = 0.0
    sumaaup = 0.0
    pre = "L"

    if doit == 1:
        pre = ""
        
        for x in range (binsaved2*m2, ca):
            a2 = 1.0*x/m2
            ff = "(%f) + %f*exp( - (x-%f)*(x-%f)/(2*(%f)*(%f)))"%(p0, a2, p2, p2, p3, p3)
            func = TF1("func", ff, 0, 50)
            sum2 = 0.0
            temporal1 = 0.0
            for i in range (1, 51):
                temporal1 = temporal1 + (h.GetBinContent(i))*log(func.Eval(i, 0, 0)) - func.Eval(i, 0, 0)*log( 2.718 ) - log( Factorial( h.GetBinContent(i)))
            if temporal1 < ( temp2 -thresholdlikelihood ) and flag == 0:
                flag = 1
                #print " Break at :",a2
                #print " Value of Chi plus: ",temporal1
                binerrup = a2
                sumaaup = temporal1

        for x in range (c0, binsaved2*m2):
            a2 = (1.0*c0)/m2 + binsaved2 - 1.0*(x)/m2
            ff = "(%f) + %f*exp( - (x-%f)*(x-%f)/(2*(%f)*(%f)))"%(p0, a2, p2, p2, p3, p3)
            func = TF1("func", ff, 0, 50)
            sum2 = 0.0
            temporal2 = 0.0     
            for i in range (1, 51):
                temporal2 = temporal2 + (h.GetBinContent(i))*log(func.Eval(i, 0, 0)) - func.Eval(i, 0, 0)*log( 2.718 ) - log( Factorial( h.GetBinContent(i)))
            if temporal2 < (temp2 -thresholdlikelihood) and flags == 0:
                flags = 1
                #print " Break at :",a2
                #print " Value of Chi minus: ",temporal2
                binerr = a2
                sumaa = temporal2
      

    ac.SetTitle("Calculating p_{1} from Likelihold")
    ac.GetYaxis().SetTitle("Likelihold")
    ac.GetXaxis().SetTitle(" p_{1} ")
    ac.SetStats(0)
    ac.SetLineColor(kBlue)
    ac.SetLineWidth(4)

    ac.Draw("")
    text = TLatex(binsaved2*0.75, temp2*0.68, " p_{1} = %.2f "%binsaved2)
    mark.DrawMarker(binsaved2, temp2)

    errm.DrawMarker(binerr, sumaa)
    errm.DrawMarker(binerrup, sumaaup)
    text.Draw("same")
    if doit == 1:
        text2 = TLatex(binsaved2*0.75, 2.9*temp2, " \Delta p_{<} = %.2f "%abs(binerr-binsaved2))
        text3 = TLatex(binsaved2*0.75, 3.3*temp2, " \Delta p_{>} = %.2f "%abs(binerrup-binsaved2))
        text2.Draw("same")
        text3.Draw("same")


    c1.SaveAs(pre+"LP1.pdf")
    

#_________________________________________________________________________
def lp0(c0 = 0, ca= 10, m2= 1.0, doit=0):
    ac =TH1F( "ac", "ac", (ca-c0), c0/m2, ca/m2)#, 500, -200, 0)

    
    temp2 = -9999.
    binsaved2 = 0
    for x in range (c0, ca):
        a2 = 1.0*x/m2
        ff = "(%f) + %f*exp( (-1.0)*(x-%f)*(x-%f)/(2*(%f)*(%f)))"%(a2, p1, p2, p2, p3, p3)
        func = TF1("func", ff, 0, 50)
        sum2 = 0.0      
        for i in range (1, 51):
            sum2 = sum2 + (h.GetBinContent(i))*log(func.Eval(i, 0, 0)) - func.Eval(i, 0, 0) - log( Factorial( h.GetBinContent(i)))
        if abs(sum2) < abs(temp2):
            temp2 = sum2
            binsaved2 = x/m2
        ac.SetBinContent(x-c0, sum2)
        #ac.Fill(a2, sum2)
    #print "L max for p0: ",temp2
    #print " value of p0: ",binsaved2
    
    binerr = 0.0
    binerrup = 0.0
    flag = 0
    flags = 0
    sumaa = 0.0
    sumaaup = 0.0
    pre = "L"

    if doit==1:
        pre = ""
        
        for x in range (binsaved2*m2, ca):
            a2 = 1.0*x/m2
            ff = "(%f) + %f*exp( (-1.0)*(x-%f)*(x-%f)/(2*(%f)*(%f)))"%(a2, p1, p2, p2, p3, p3)
            func = TF1("func", ff, 0, 50)
            sum2 = 0.0
            temporal1 = 0.0
            for i in range (1, 51):
                temporal1 = temporal1 + (h.GetBinContent(i))*log(func.Eval(i, 0, 0)) - func.Eval(i, 0, 0)*log( 2.718 ) - log( Factorial( h.GetBinContent(i)))
            if temporal1 < ( temp2 -thresholdlikelihood ) and flag == 0:
                flag = 1
                #print " Break at :",a2
                #print " Value of Chi plus: ",temporal1
                binerrup = a2
                sumaaup = temporal1

        for x in range (c0, binsaved2*m2):
            a2 = (1.0*c0)/m2 + binsaved2 - 1.0*(x)/m2
            ff = "(%f) + %f*exp( (-1.0)*(x-%f)*(x-%f)/(2*(%f)*(%f)))"%(a2, p1, p2, p2, p3, p3)
            func = TF1("func", ff, 0, 50)
            sum2 = 0.0
            temporal2 = 0.0     
            for i in range (1, 51):
                temporal2 = temporal2 + (h.GetBinContent(i))*log(func.Eval(i, 0, 0)) - func.Eval(i, 0, 0)*log( 2.718 ) - log( Factorial( h.GetBinContent(i)))
            if temporal2 < (temp2 -thresholdlikelihood) and flags == 0:
                flags = 1
                #print " Break at :",a2
                #print " Value of Chi minus: ",temporal2
                binerr = a2
                sumaa = temporal2
      

    ac.SetTitle("Calculating p_{0} from Likelihold")
    ac.GetYaxis().SetTitle("Likelihold")
    ac.GetXaxis().SetTitle(" p_{0} ")
    ac.SetStats(0)
    ac.SetLineColor(kRed)
    ac.SetLineWidth(4)
    ac.Draw("")
    
    text = TLatex(binsaved2*2.75, 1*temp2, " p_{0} = %.2f "%binsaved2)
    mark.DrawMarker(binsaved2, temp2)

    errm.DrawMarker(binerr, sumaa)
    errm.DrawMarker(binerrup, sumaaup)

    text.Draw("same")

    if doit == 1:
        text2 = TLatex(binsaved2*0.75, 4*temp2, " \Delta p_{<} = %.2f "%abs(binerr-binsaved2))
        text3 = TLatex(binsaved2*0.75, 4.8*temp2, " \Delta p_{>} = %.2f "%abs(binerrup-binsaved2))
        text2.Draw("same")
        text3.Draw("same")


    c1.SaveAs(pre+"LP0.pdf")



#_________________________________________________________________________
def lp2p3(c0 = 0, ca= 10, m2= 1.0, a0= 0, a1 = 10, m3 = 1.0, doit=0):
    ak =TH2F( "ak", "ak", ca-c0, c0/m2, ca/m2, a1-a0, a0/m3, a1/m3)
    temp2 = 9999.
    binsaved2 = 0.0
    binsaved3 = 0.0
    for x2 in range (c0, ca):
        b0 = 1.0*x2/m2
        for y in range (a0, a1):
            a2 = 1.0*y/m3
            ff = "(%f) + %f*exp( (-1.0)*(x-%f)*(x-%f)/(2*(%f)*(%f)))"%(p0, p1, b0, b0, a2, a2) 
            func = TF1("func", ff, 0, 41)
            sum2 = 0.0
            for i in range (1, 41):
                sum2 = sum2 + (h.GetBinContent(i))*log(func.Eval(i, 0, 0)) - func.Eval(i, 0, 0) - log( Factorial( h.GetBinContent(i)))
            if abs(sum2) < abs(temp2):
                temp2 = sum2
                binsaved2 = b0
                binsaved3 = a2

            ak.SetBinContent(x2 - c0, y-a0, sum2) 

    #print "max Likelihood for p2p3: ",temp2
    #print " value of p2: ",binsaved2
    #print " Value of p3: ",binsaved3

    ak.SetTitle("Calculating p_{3} and p_{2} from Likelihood")
    ak.GetYaxis().SetTitle(" p_{3}")
    ak.GetXaxis().SetTitle("  p_{2} ")
    ak.SetStats(0)

    ak.Draw("colz")   
    mark.DrawMarker(binsaved2, binsaved3)
    
    pre = "L"
    scale1= 1.3
    scale2 = 0.7
    
    if doit==1:
        scale1= 0.9
        scale2=0.88
        pre = ""
        for x2 in range (c0, ca):
            b0 = 1.0*x2/m2
            for y in range (a0, a1):
                a2 = 1.0*y/m3
                ff ="(%f) + %f*exp( (-1.0)*(x-%f)*(x-%f)/(2*(%f)*(%f)))"%(p0, p1, b0, b0, a2, a2) 
                func = TF1("func", ff, 0, 41)
                sum2 = 0.0
                for i in range (1, 41):
                    sum2 = sum2 + (h.GetBinContent(i))*log(func.Eval(i, 0, 0)) - func.Eval(i, 0, 0) - log( Factorial( h.GetBinContent(i)))
                if abs(sum2) <= abs(temp2 - thresholdlikelihood):
                    if abs(sum2) >= abs(temp2 - 0.93*thresholdlikelihood):
                        mk.DrawMarker(b0, a2)
                    

    text2 = TLatex(0.96*binsaved2, binsaved3*scale1, " P2: %.2f"%binsaved2)
    text3 = TLatex(0.96*binsaved2, binsaved3*scale2, " P3: %.2f"%binsaved3)    
    text2.Draw("same")
    text3.Draw("same")  
    c1.SaveAs(pre+"LP2P3.pdf")


#_________________________________________________________________________
def lp1p2(c0 = 0, ca= 10, m2= 1.0, a0= 0, a1 = 10, m3 = 1.0, doit=0):
    ak =TH2F( "ak", "ak", ca-c0, c0/m2, ca/m2, a1-a0, a0/m3, a1/m3)
    temp2 = 9999.
    binsaved2 = 0.0
    binsaved3 = 0.0
    for x2 in range (c0, ca):
        b0 = 1.0*x2/m2 
        for y in range (a0, a1):
            a2 = 1.0*y/m3
            ff = "(%f) + %f*exp( (-1)*(x-%f)*(x-%f)/(2*(%f)*(%f)))"%(p0, b0, a2, a2, p3, p3) 
            func = TF1("func", ff, 0, 41)
            sum2 = 0.0
            for i in range (1, 41):
                sum2 = sum2 + (h.GetBinContent(i))*log(func.Eval(i, 0, 0)) - func.Eval(i, 0, 0) - log( Factorial( h.GetBinContent(i)))
            if abs(sum2) < abs(temp2):
                temp2 = sum2
                binsaved2 = b0
                binsaved3 = a2

            ak.SetBinContent(x2 - c0, y-a0, sum2) 

    #print "max Likelihood for p2: ",temp2
    #print " value of p1: ",binsaved2
    #print " Value of p2: ",binsaved3 
    ak.SetTitle("Calculating p_{1} and p_{2} from Likelihood")
    ak.GetYaxis().SetTitle(" p_{2}")
    ak.GetXaxis().SetTitle("  p_{1} ")
    ak.SetStats(0)
    ak.SetMaximum(-90)

    ak.Draw("colz")   
    mark.DrawMarker(binsaved2, binsaved3)

    pre = "L"
    scale1= 1.3
    scale2 = 0.7
    
    if doit==1:
        scale1= 0.93
        scale2=0.91
        pre = ""
        for x2 in range (c0, ca):
            b0 = 1.0*x2/m2
            for y in range (a0, a1):
                a2 = 1.0*y/m3
                ff ="(%f) + %f*exp( (-1)*(x-%f)*(x-%f)/(2*(%f)*(%f)))"%(p0, b0, a2, a2, p3, p3) 
                func = TF1("func", ff, 0, 41)
                sum2 = 0.0
                for i in range (1, 41):
                    sum2 = sum2 + (h.GetBinContent(i))*log(func.Eval(i, 0, 0)) - func.Eval(i, 0, 0) - log( Factorial( h.GetBinContent(i)))
                if abs(sum2) <= abs(temp2 - thresholdlikelihood):
                    if abs(sum2) >= abs(temp2 - 0.93*thresholdlikelihood):
                        mk.DrawMarker(b0, a2)      

    text2 = TLatex(0.96*binsaved2, binsaved3*scale1, " P1: %.2f"%binsaved2)
    text3 = TLatex(0.96*binsaved2, binsaved3*scale2, " P2: %.2f"%binsaved3)    
    text2.Draw("same")
    text3.Draw("same")           

    c1.SaveAs(pre+"LP1P2.pdf")

#_________________________________________________________________________________________________________________

def lp1p3(c0 = 0, ca= 10, m2= 1.0, a0= 0, a1 = 10, m3 = 1.0, doit=0):
    ak =TH2F( "ak", "ak", ca-c0, c0/m2, ca/m2, a1-a0, a0/m3, a1/m3)
    temp2 = 9999.
    binsaved2 = 0.0
    binsaved3 = 0.0
    for x2 in range (c0, ca):
        b0 = 1.0*x2/m2
        for y in range (a0, a1):
            a2 = 1.0*y/m3
            ff = "(%f) + %f*exp( (-1)*(x-%f)*(x-%f)/(2*(%f)*(%f)))"%(p0, b0, p2, p2, a2, a2) 
            func = TF1("func", ff, 0, 41)
            sum2 = 0.0
            for i in range (1, 41):
                sum2 = sum2 +(h.GetBinContent(i))*log(func.Eval(i, 0, 0)) - func.Eval(i, 0, 0) - log( Factorial( h.GetBinContent(i)))
            if abs(sum2) < abs(temp2):

                temp2 = sum2
                binsaved2 = b0
                binsaved3 = a2
                
            ak.SetBinContent(x2 - c0, y-a0, sum2) 


    #print "max Likelihood for p2: ",temp2
    #print " value of p1: ",binsaved2
    #print " Value of p3: ",binsaved3
    ak.SetTitle("Calculating p_{1} and p_{3} from Likelihood")
    ak.GetYaxis().SetTitle(" p_{3}")
    ak.GetXaxis().SetTitle("  p_{1} ")
    ak.SetStats(0)

    ak.Draw("colz")   
    mark.DrawMarker(binsaved2, binsaved3)

    pre = "L"
    scale1= 1.3
    scale2 = 0.7
    
    if doit==1:
        scale1= 0.91
        scale2=0.88
        pre = ""
        for x2 in range (c0, ca):
            b0 = 1.0*x2/m2
            for y in range (a0, a1):
                a2 = 1.0*y/m3
                ff ="(%f) + %f*exp( (-1)*(x-%f)*(x-%f)/(2*(%f)*(%f)))"%(p0, b0, p2, p2, a2, a2) 
                func = TF1("func", ff, 0, 41)
                sum2 = 0.0
                for i in range (1, 41):
                    sum2 = sum2 + (h.GetBinContent(i))*log(func.Eval(i, 0, 0)) - func.Eval(i, 0, 0) - log( Factorial( h.GetBinContent(i)))
                if abs(sum2) <= abs(temp2 - thresholdlikelihood):
                    if abs(sum2) >= abs(temp2 - 0.93*thresholdlikelihood):
                        mk.DrawMarker(b0, a2)
                    
    text2 = TLatex(0.96*binsaved2, binsaved3*scale1, " P1: %.2f"%binsaved2)
    text3 = TLatex(0.96*binsaved2, binsaved3*scale2, " P3: %.2f"%binsaved3)    
    text2.Draw("same")
    text3.Draw("same")    

    c1.SaveAs(pre+"LP1P3.pdf")


#_________________________________________________________________________________________________________________

def lp0p1(c0 = 0, ca= 10, m2= 1.0, a0= 0, a1 = 10, m3 = 1.0, doit=0):

    ak =TH2F( "ak", "ak", ca-c0, c0/m2, ca/m2, a1-a0, a0/m3, a1/m3)
    temp2 = 9999.
    binsaved2 = 0.0
    binsaved3 = 0.0
    
    for x2 in range (c0, ca):
        b0 = 1.0*x2/m2
        for y in range (a0, a1):
            a2 = 1.0*y/m3
            ff = "(%f) + %f*exp( (-1)*(x-%f)*(x-%f)/(2*(%f)*(%f)))"%(b0, a2, p2, p2, p3, p3) 

            func = TF1("func", ff, 0, 41)
            sum2 = 0.0
            
            for i in range (1, 41):
                sum2 = sum2 +(h.GetBinContent(i))*log(func.Eval(i, 0, 0)) - func.Eval(i, 0, 0) - log( Factorial( h.GetBinContent(i)))
            if abs(sum2) < abs(temp2):

                temp2 = sum2
                binsaved2 = b0
                binsaved3 = a2
                
            ak.SetBinContent(x2 - c0, y-a0, sum2) 

    #print "max Likelihood for p0p1: ",temp2
    #print " value of p0: ",binsaved2
    #print " Value of p1: ",binsaved3
 
    ak.SetTitle("Calculating p_{0} and p_{1} from Likelihood")
    ak.GetYaxis().SetTitle(" p_{1}")
    ak.GetXaxis().SetTitle("  p_{0} ")
    ak.SetStats(0)
    ak.Draw("colz")   
    mark.DrawMarker(binsaved2, binsaved3)
    
    pre = "L"
    scale1= 1.3
    scale2 = 0.7
    
    if doit==1:
        pre = ""
        scale1= 1.02
        scale2= 1
        for x2 in range (c0, ca):
            b0 = 1.0*x2/m2             
            for y in range (a0, a1):
                a2 = 1.0*y/m3
                ff ="(%f) + %f*exp( (-1)*(x-%f)*(x-%f)/(2*(%f)*(%f)))"%(b0, a2, p2, p2, p3, p3) 
                func = TF1("func", ff, 0, 41)
                sum2 = 0.0        
                for i in range (1, 41):
                    sum2 = sum2 + (h.GetBinContent(i))*log(func.Eval(i, 0, 0)) - func.Eval(i, 0, 0) - log( Factorial( h.GetBinContent(i)))
                if abs(sum2) <= abs(temp2 - thresholdlikelihood):
                    if abs(sum2) >= abs(temp2 - 0.93*thresholdlikelihood):
                        mk.DrawMarker(b0, a2)
        
    text2 = TLatex(0.56*binsaved2, binsaved3*scale1, " P0^{2}: %.2f"%binsaved2)
    text3 = TLatex(0.56*binsaved2, binsaved3*scale2, " P1: %.2f"%binsaved3)    
    text2.Draw("same")
    text3.Draw("same")  
                    
    c1.SaveAs(pre+"LP0P1.pdf")


    

#######################################################################


#lp0p1(100, 300, 100, 240, 290, 10, 1)
#lp0p1(0, 100, 10, 0, 100, 1, 0)
#lp1p3(230, 300, 10, 300, 600, 100, 1)
#lp1p3(0, 60, 1, 0, 120, 10, 0)
#lp1p2(0, 100, 1, 0, 60, 1, 0)
#lp1p2(2100, 3000, 100, 1900, 2300, 100, 1)
#lp2p3(1900, 2300, 100.,400,600, 100, 1)   
#lp2p3(0, 80, 1 ,0 , 120, 10, 0)
'''
lp0(0, 3600, 1000.,1)
lp1(2000, 3500, 100.,1)
lp2(1700, 2500, 100.,1)
lp3(4000, 6000, 1000.,1)
'''

#lp0(0, 3000, 100.,1)
#lp1(0, 7000, 100.,1)
#lp2(0, 7000, 100.,1)
lp3(0, 700000, 1.,1)

#######################################################################


#dop0p1(100, 400, 100, 240, 300, 10, 1)
#dop0p1(10, 100, 10, 10, 800, 10, 0)
#dop1p3(220, 300, 10, 400, 600, 100,1)
#dop1p3(0, 600, 10, 0, 100, 10,0)
#dop1p2(2300, 3100, 100, 2000, 2200, 100, 1)
#dop1p2(0, 800, 10, 0, 500, 10, 0)
#dop2p3(2000, 2300, 100.,400,600, 100, 1)
#dop2p3(0, 600, 10.,0,120, 10, 0)
#dop0(500, 3600, 1000., 0)
#dop0(0, 3000, 100., 1)

#dop1(2000, 3500, 100., 1)
#dop1(500, 7000, 100., 1)
#dop2(1700, 2500, 100.)
#dop2(0, 7000, 100., 1)
#dop3(4000, 6000, 1000., 0)
#dop3(0, 7000, 100., 0)
#minimize("p3")
