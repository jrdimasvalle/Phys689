import sys
from ROOT import *


c1 = TCanvas()
c1.SetGridx()
c1.SetGridy()
c1.SetTickx()
c1.SetTicky()


f = TFile("lab3.root")
h = f.Get("h1")
h2 = f.Get("h2")
h3 = f.Get("h3")

bins_number = h.GetSize()-2
max_bin_value = 5.0
alpha = 0.5
mass_0 = 3.1
sigma = 0.3


marker_position = TMarker(1, 40, 34)

marker_1sigma = TMarker(1, 40, 33)
marker_2sigma = TMarker(1, 40, 24)

marker_cl_2sigma = TMarker(5, 40, 7)
marker_cl_2sigma.SetMarkerColor(2)
marker_cl_1sigma = TMarker(6, 40, 5)
marker_cl_1sigma.SetMarkerColor(2)
thresholdchisquared = 1
thresholdlikelihood = 1


#___________________________________________________
def Factorial( n = 1 ):
    if (n < 0): return 1
    x= 1.0
    for b in range (1, n+1):
        x = x*b
    return x

        
#_________________________________________________________________________
def ChiSquare_TwoVariables(min_limit_first = 0, max_limit_first= 10, divisions_per_bin = 1.0, min_limit_second = 0 , max_limit_second = 10, divisions_per_bin_second = 1.0, doit = 0):
    
    ChiSquare_Histogram =TH2F("ChiSquare_Histogram", "ChiSquare_Histogram", max_limit_first - min_limit_first, min_limit_first/divisions_per_bin, max_limit_first/divisions_per_bin, max_limit_second - min_limit_second, min_limit_second/divisions_per_bin_second, max_limit_second/divisions_per_bin_second)

    f1 = TH1F("f1", "f1", 50, 0, 5)
    chi_square = 99999.
    position_x = 0.0
    position_y = 0.0


    for x1 in range (min_limit_first, max_limit_first+1):
        xx = 1.0*x1/divisions_per_bin

 
        for y in range (min_limit_second, max_limit_second+1):
            yy = 1.0*y/divisions_per_bin

   
            ff = "(%f)*exp((-1.0)*(%f)*(x)) + (%f)*exp((-1.0)*(x - %f)*(x-%f)/(%f*%f*2))"%(xx, alpha, yy, mass_0, mass_0, sigma, sigma)
            function = TF1("function", ff, 0, bins_number)
            chi_square_temp = 0.0
          
            for ii in range (1, bins_number+1):
                f_value = 1.0*function.Eval(1.0*max_bin_value*ii/bins_number, 0, 0)
                chi_square_temp = chi_square_temp + (1.0*f_value - 1.0*h.GetBinContent(ii))*(1.0*f_value - 1.0*h.GetBinContent(ii))/(1.0*f_value)

            ChiSquare_Histogram.SetBinContent(x1 - min_limit_first, y-min_limit_second, chi_square_temp)

             
            if chi_square_temp < chi_square:
                
                chi_square = chi_square_temp
                position_x = xx
                position_y = yy

        
    print "Value of ChiSquare: ",chi_square
    print "Value of first variable: ",position_x
    print "Value of second variable: ",position_y

    
    

    ChiSquare_Histogram.SetTitle("Calculating B and S from \Chi^{2}")
    ChiSquare_Histogram.GetYaxis().SetTitle("S [Data]")
    ChiSquare_Histogram.GetXaxis().SetTitle("B [Background]")
    ChiSquare_Histogram.SetStats(0)

    ChiSquare_Histogram.Draw("colz")   
    marker_position.DrawMarker(position_x, position_y)
 


    if doit==1:
        pre = ""
        for x1 in range (min_limit_first, max_limit_first+1):
            xx = 1.0*x1/divisions_per_bin

     
            for y in range (min_limit_second, max_limit_second+1):
                yy = 1.0*y/divisions_per_bin

       
                ff = "(%f)*exp((-1.0)*(%f)*(x)) + (%f)*exp((-1.0)*(x - %f)*(x-%f)/(%f*%f*2))"%(xx, alpha, yy, mass_0, mass_0, sigma, sigma)
                function = TF1("function", ff, 0, bins_number)
                chi_square_temp = 0.0
              
                for ii in range (1, bins_number+1):
                    f_value = 1.0*function.Eval(1.0*max_bin_value*ii/bins_number, 0, 0)
                    chi_square_temp = chi_square_temp + (1.0*f_value - 1.0*h.GetBinContent(ii))*(1.0*f_value - 1.0*h.GetBinContent(ii))/(1.0*f_value)

                if abs(chi_square_temp) < abs(chi_square)+1:
                    if abs(chi_square_temp) > (abs(chi_square)) + 0.89:
                        marker_1sigma.DrawMarker(xx, yy)


                if abs(chi_square_temp) < abs(chi_square)+2:
                    if abs(chi_square_temp) > (abs(chi_square)) + 1.89:
                        marker_2sigma.DrawMarker(xx, yy)


    text2 = TLatex(0.96*position_x, position_y*1.025, " S: %.2f"%position_y)
    text3 = TLatex(0.96*position_x, position_y*1.012, " B: %.2f"%position_x)    
    text2.Draw("same")
    text3.Draw("same") 
    c1.SaveAs(pre+"Chi_Square.pdf")


#_________________________________________________________________________________
def Likelihood_second(A_0, hx = TH1F("ac", "ac", 10, 0, 10), is_constant = 0):

    nmax = hx.GetSize()-1
    sum2 = 0.0      
    for i in range (1, nmax):


        ff = "(%f)*exp((-1.0)*(%f)*(x)) + (%f)*exp((-1.0)*(x - %f)*(x-%f)/(%f*%f*2))"%(42.8, alpha, 31.1, mass_0, mass_0, sigma, sigma)
        function = TF1("function", ff, 0, 41)

        f_value = 1.0*function.Eval(1.0*max_bin_value*i/bins_number, 0, 0)

        if is_constant == 1:
            f_value = A_0        
        
        sum2 = sum2 + (hx.GetBinContent(i))*log(f_value) - f_value - log( Factorial( hx.GetBinContent(i)))

        

    return sum2




#_________________________________________________________________________
def Likelihood_TwoVariables(min_limit_first = 0, max_limit_first= 10, divisions_per_bin = 1.0, min_limit_second = 0 , max_limit_second = 10, divisions_per_bin_second = 1.0, doit = 0):
    
    Likelihood_Histogram =TH2D("Likelihood_Histogram", "Likelihood_Histogram", max_limit_first - min_limit_first, min_limit_first/divisions_per_bin, max_limit_first/divisions_per_bin, max_limit_second - min_limit_second, min_limit_second/divisions_per_bin_second, max_limit_second/divisions_per_bin_second)

   
    likelihood = -99999.
    position_x = 0.0
    position_y = 0.0


    normalization = 0.0
                            
    for x1 in range (min_limit_first, max_limit_first+1):
        xx = 1.0*x1/divisions_per_bin
        for y in range (min_limit_second, max_limit_second+1):
            yy = 1.0*y/divisions_per_bin
            ff = "(%f)*exp((-1.0)*(%f)*(x)) + (%f)*exp((-1.0)*(x - %f)*(x-%f)/(%f*%f*2))"%(xx, alpha, yy, mass_0, mass_0, sigma, sigma)
            function = TF1("function", ff, 0, 41)
            likelihood_temp = 0.0
          
            for ii in range (1, bins_number+1):
                f_value = 1.0*function.Eval(1.0*max_bin_value*ii/bins_number, 0, 0)
                likelihood_temp = likelihood_temp + (h.GetBinContent(ii))*log(f_value) - f_value - log( Factorial( h.GetBinContent(ii)))


            normalization = normalization + abs(likelihood_temp)
            Likelihood_Histogram.SetBinContent(x1 - min_limit_first, y-min_limit_second, likelihood_temp)

            if abs(likelihood_temp) < abs(likelihood):
                
                likelihood = likelihood_temp
                position_x = xx
                position_y = yy




    print "Normalization: ",normalization
    cl_value = 0.95
    looked_value = cl_value*normalization
    tempo = 0.0
    for y in range (0,max_limit_first+1):
        for x1 in range(0,2*max_limit_second+1):
            tempo = tempo +  abs(Likelihood_Histogram.GetBinContent(Likelihood_Histogram.GetBin(x1, y)))
            
        if abs(tempo) > abs(looked_value):
            print "CL 95 found at y: ",1.0*(1.0*y+1.0*min_limit_second)/divisions_per_bin
            break

    print "Value of Likelihood: ",likelihood
    print "Value of first variable: ",position_x
    print "Value of second variable: ",position_y

    Likelihood_Histogram.Draw("colz")
    
    pre = "Long"
    Likelihood_Histogram.SetTitle("Calculating B and S from log [ Likelihoodd ]")
    Likelihood_Histogram.GetYaxis().SetTitle("S [Data]")
    Likelihood_Histogram.GetXaxis().SetTitle("B [Background]")
    Likelihood_Histogram.SetStats(0)

    Likelihood_Histogram.Draw("colz")   
    marker_position.DrawMarker(position_x, position_y)


    

    if (doit == 1):
        pre = ""
        for x1 in range (min_limit_first, max_limit_first+1):
            xx = 1.0*x1/divisions_per_bin
            for y in range (min_limit_second, max_limit_second+1):
                yy = 1.0*y/divisions_per_bin
                ff = "(%f)*exp((-1.0)*(%f)*(x)) + (%f)*exp((-1.0)*(x - %f)*(x-%f)/(%f*%f*2))"%(xx, alpha, yy, mass_0, mass_0, sigma, sigma)
                function = TF1("function", ff, 0, bins_number)
                likelihood_temp = 0.0
              
                for ii in range (1, bins_number+1):
                    f_value = 1.0*function.Eval(1.0*max_bin_value*ii/bins_number, 0, 0)
                    likelihood_temp = likelihood_temp + (h.GetBinContent(ii))*log(f_value) - f_value - log( Factorial( h.GetBinContent(ii)))

       
                        
                if abs(likelihood_temp)< abs(likelihood)+thresholdlikelihood:
                    if abs(likelihood_temp) > abs(likelihood)+thresholdlikelihood - 0.1:
                        marker_1sigma.DrawMarker(xx,yy)
                        

                if abs(likelihood_temp)< abs(likelihood)+thresholdlikelihood + 1:
                    if abs(likelihood_temp) > abs(likelihood)+thresholdlikelihood +1 - 0.1:
                        marker_2sigma.DrawMarker(xx,yy)
                        

    text2 = TLatex(0.96*position_x, position_y*1.06, " S: %.2f"%position_y)
    text3 = TLatex(0.96*position_x, position_y*1.034, " B: %.2f"%position_x)    
    text2.Draw("same")
    text3.Draw("same") 


    c1.SaveAs(pre+"aaLikelihood.pdf")





#_____________________________________________________________________________________________________
def Pseudo_Experiment(number_of_experiments = 10, hx = TH1D("fakeasdfas", "fakeasfdasa", 50, 0, 5), is_constant = 1):


    A_0= Likelihood_OneVariable(hx, 0, 1000, 10, 0, is_constant)
    position_x = Likelihood_second(A_0, hx, is_constant)

    print A_0
    
    fake_experiment = TH1D("fake", "fake", 10, 0, 10)
    Likelihood_Histogram = TH1D("ll", "ll", 200, int(position_x*1.5), int(position_x*0.5))




    for j in range (0, number_of_experiments):

        binn = 10
        if is_constant == 0:
            binn = 50
    
        for k in range(0, binn):
            
            if is_constant == 1:
                 xx = TRandom1().PoissonD(A_0)
                 fake_experiment.SetBinContent(k+1, xx)
                 
            if is_constant == 0:
                 xx = TRandom1().PoissonD(mass_0)
                 ff = "(42.80)*exp((-1.0)*(%f)*(x)) + (31.8)*exp((-1.0)*(x - %f)*(x-%f)/(%f*%f*2))"%(alpha, mass_0, mass_0, sigma, sigma)
                 function = TF1("func", ff, 0, 50)


                 fake_experiment.SetBinContent(k+1, function.Eval(xx, 0,0))


        


        if is_constant == 1:
                Likelihood_Histogram.Fill(Likelihood_second(A_0, fake_experiment))

        if is_constant == 0:

                Likelihood_Histogram.Fill(Likelihood_second(31.8, fake_experiment))            

    
    name =  hx.GetTitle()
    Likelihood_Histogram.SetTitle("Likelihood from PseudoExperiments on "+name)
    Likelihood_Histogram.GetYaxis().SetTitle(" Frequency")
    Likelihood_Histogram.GetXaxis().SetTitle(" log [Likelihood] ")
    Likelihood_Histogram.SetStats(0)
    Likelihood_Histogram.SetLineColor(kOrange)
    Likelihood_Histogram.SetLineWidth(4)
    Likelihood_Histogram.Draw("")


        
    Likelihood_Histogram.Draw()
    

  
    bin_x = Likelihood_Histogram.GetXaxis().FindBin(position_x)
    maximum_y = Likelihood_Histogram.GetMaximum()

    for j in range(0, 10*maximum_y):      
        marker_cl_2sigma.DrawMarker(position_x , 0.1*j)


    #total_integral = Likelihood_Histogram.Integral(1, 99999)


    #partial_integral = Likelihood_Histogram.Integral(1, bin_x)/total_integral



    #text3 = TLatex(1.6*position_x, 0.5*maximum_y, "p value: %.3f"%partial_integral)    

    #text3.Draw("same")
    #print "A_0: ",A_0
    
    c1.SaveAs("fake"+name+".pdf")

    
#____________________________________________-
def Likelihood_OneVariable(hx = TH1D("fakeasdfas", "fakeasfdasa", 50, 0, 50), lower_limit= 0, upper_limit = 10, divisions_per_bin= 1.0, doit=0, is_constant = 1):

    Likelihood_Histogram =TH1D( "Likelihood_Histogram", "Likelihood_Histogram", (upper_limit - lower_limit), lower_limit/divisions_per_bin, upper_limit/divisions_per_bin)#, 500, -200, 0)
    real_Likelihood_Histogram =TH1D( "real_Likelihood_Histogram", "real_Likelihood_Histogram", (upper_limit - lower_limit), lower_limit/divisions_per_bin, upper_limit/divisions_per_bin)#, 500, -200, 0)
    likelihood = - 99999.
    position_x = 0

    bins_number = hx.GetSize()-2
    
    for x1 in range (lower_limit, upper_limit+1):
        xx = 1.0*x1/divisions_per_bin
        ff = "(42.80)*exp((-1.0)*(%f)*(x)) + (%f)*exp((-1.0)*(x - %f)*(x-%f)/(%f*%f*2))"%(alpha, xx, mass_0, mass_0, sigma, sigma)
        function = TF1("func", ff, 0, bins_number)
        likelihood_temp = 0.0
        for ii in range (1, bins_number+1):
                f_value = 1.0*function.Eval(1.0*max_bin_value*ii/bins_number, 0, 0)

                if is_constant == 1:
                    f_value = xx
                likelihood_temp = likelihood_temp + (hx.GetBinContent(ii))*log(f_value) - f_value - log( Factorial( hx.GetBinContent(ii)))

        Likelihood_Histogram.SetBinContent(x1 - lower_limit, likelihood_temp)
        real_Likelihood_Histogram.SetBinContent(x1 - lower_limit, 1.0*exp(likelihood_temp))
        if abs(likelihood_temp) < abs(likelihood):    
            likelihood = likelihood_temp
            position_x = xx

 
    print likelihood
    
    bin_err_up = 0.0
    bin_err_up_2sigma = 0.0
    bin_err_dw_2sigma = 0.0
    bin_err_dw = 0.0


    has_it_up = 0
    has_it_dw = 0

    has_it_up_2sigma = 0
    has_it_dw_2sigma = 0

    likelihood_up = 0.0
    likelihood_up_2sigma = 0.0
    likelihood_dw = 0.0
    likelihood_dw_2sigma = 0.0
    pre = "Long"

    name = hx.GetTitle()

    if doit==1:
        pre = ""


        for x1 in range (position_x*divisions_per_bin, upper_limit):
            xx = 1.0*x1/divisions_per_bin
            ff = "(42.80)*exp((-1.0)*(%f)*(x)) + (%f)*exp((-1.0)*(x - %f)*(x-%f)/(%f*%f*2))"%(alpha, xx, mass_0, mass_0, sigma, sigma)
            function = TF1("func", ff, 0, bins_number)
            likelihood_temp = 0.0
            for ii in range (1, bins_number+1):
                    f_value = 1.0*function.Eval(1.0*max_bin_value*ii/bins_number, 0, 0)

                    if is_constant == 1:
                        f_value = xx
                    likelihood_temp = likelihood_temp + (hx.GetBinContent(ii))*log(f_value) - f_value - log( Factorial( hx.GetBinContent(ii)))
                    
            if abs(likelihood_temp)> abs(likelihood)+thresholdlikelihood and has_it_up == 0:
                has_it_up = 1
                bin_err_up = xx
                likelihood_up = likelihood_temp

            if abs(likelihood_temp)> abs(likelihood)+thresholdlikelihood +1  and has_it_up_2sigma == 0:
                has_it_up_2sigma = 1
                bin_err_up_2sigma = xx
                likelihood_up_2sigma = likelihood_temp     
               
            
        for x1 in range (lower_limit, position_x*divisions_per_bin):
            xx = 1.0*x1/divisions_per_bin
            ff = "(42.80)*exp((-1.0)*(%f)*(x)) + (%f)*exp((-1.0)*(x - %f)*(x-%f)/(%f*%f*2))"%(alpha, xx, mass_0, mass_0, sigma, sigma)
            function = TF1("func", ff, 0, bins_number)
            likelihood_temp = 0.0
            for ii in range (1, bins_number+1):
                    f_value = 1.0*function.Eval(1.0*max_bin_value*ii/bins_number, 0, 0)
                    if is_constant == 1:
                        f_value = xx
                    likelihood_temp = likelihood_temp + (hx.GetBinContent(ii))*log(f_value) - f_value - log( Factorial( hx.GetBinContent(ii)))

            if abs(likelihood_temp)< abs(likelihood)+thresholdlikelihood and has_it_dw == 0:
                has_it_dw = 1
                bin_err_dw = xx
                likelihood_dw = likelihood_temp
                
            if abs(likelihood_temp)< abs(likelihood)+thresholdlikelihood + 1and has_it_dw_2sigma == 0:
                has_it_dw_2sigma = 1
                bin_err_dw_2sigma = xx
                likelihood_dw_2sigma = likelihood_temp
                    
        print "Likelihood: ",likelihood
        print "X position: ",position_x
        print "Error up x: ",bin_err_up
        print "Likelihood err up: ",likelihood_up
        print "Error dw x: ",bin_err_dw
        print "Likelihood err dw: ",likelihood_dw

        Likelihood_Histogram.SetTitle("Calculating S from Likelihold on "+name)
        Likelihood_Histogram.GetYaxis().SetTitle("log [ Likelihold]")
        Likelihood_Histogram.GetXaxis().SetTitle(" S [Data] ")
        Likelihood_Histogram.SetStats(0)
        Likelihood_Histogram.SetLineColor(kOrange)
        Likelihood_Histogram.SetLineWidth(4)
        Likelihood_Histogram.Draw("")
        text = TLatex(position_x, likelihood, " S = %.2f "%position_x)
        
        marker_position.DrawMarker(position_x, likelihood)
        text.Draw("same")


        marker_1sigma.DrawMarker(bin_err_up, likelihood_up)
        marker_1sigma.DrawMarker(bin_err_dw, likelihood_dw)

        marker_2sigma.DrawMarker(bin_err_up_2sigma, likelihood_up_2sigma)
        marker_2sigma.DrawMarker(bin_err_dw_2sigma, likelihood_dw_2sigma)

        
        c1.SaveAs(pre+"Likelihood_1var"+name+".pdf")


        
        real_Likelihood_Histogram.SetTitle("Calculating S from Likelihold on "+name)
        real_Likelihood_Histogram.GetYaxis().SetTitle("Likelihold")
        real_Likelihood_Histogram.GetXaxis().SetTitle(" S [Data] ")
        real_Likelihood_Histogram.SetStats(0)
        real_Likelihood_Histogram.SetLineColor(kOrange)
        real_Likelihood_Histogram.SetLineWidth(4)
        real_Likelihood_Histogram.Draw("")
        text2 = TLatex(position_x, exp(likelihood), " S = %.2f "%position_x)
        marker_position.DrawMarker(position_x, exp(likelihood))
        marker_1sigma.DrawMarker(bin_err_up, exp(likelihood_up))
        marker_1sigma.DrawMarker(bin_err_dw, exp(likelihood_dw))
        marker_2sigma.DrawMarker(bin_err_up_2sigma, exp(likelihood_up_2sigma))
        marker_2sigma.DrawMarker(bin_err_dw_2sigma, exp(likelihood_dw_2sigma))
     
        text2.Draw("same")

        c1.SaveAs(pre+"real_Likelihood_1var"+name+".pdf")

       

        total_integral = real_Likelihood_Histogram.Integral(0, 999999999)
        real_Likelihood_Histogram.Scale(1/total_integral)
        
        real_Likelihood_Histogram.SetTitle("Calculating S from Normalized Likelihood on "+name)
        real_Likelihood_Histogram.GetYaxis().SetTitle("Likelihold (Normalized) ")

        axis = TGaxis()
        axis.SetMaxDigits(2)
        
        real_Likelihood_Histogram.GetXaxis().SetTitle(" S [Data] ")
        real_Likelihood_Histogram.SetStats(0)
        real_Likelihood_Histogram.SetLineColor(kOrange)
        real_Likelihood_Histogram.SetLineWidth(4)
        real_Likelihood_Histogram.Draw("")
        text2 = TLatex(position_x, exp(likelihood)*(1/total_integral), " S = %.2f "%position_x)


        marker_position.DrawMarker(position_x, exp(likelihood)*(1/total_integral))
        marker_1sigma.DrawMarker(bin_err_up, exp(likelihood_up)*(1/total_integral))
        marker_1sigma.DrawMarker(bin_err_dw, exp(likelihood_dw)*(1/total_integral))
        marker_2sigma.DrawMarker(bin_err_up_2sigma, exp(likelihood_up_2sigma)*(1/total_integral))
        marker_2sigma.DrawMarker(bin_err_dw_2sigma, exp(likelihood_dw_2sigma)*(1/total_integral))
        text2.Draw("same")


        total_integral = real_Likelihood_Histogram.Integral(0, 999999999)
        print "Total: ",total_integral
        central_value = int(position_x*divisions_per_bin)
        flag_1 = 0
        flag_2 = 0
        CL_dw = 0.0
        CL_up = 0.0
        
        for k in range(lower_limit, central_value):
            value_int = real_Likelihood_Histogram.Integral(1, k - lower_limit + 1)
            if value_int > 0.025 and flag_1 == 0:    
                print " 95% confidence lower limit: ",1.0*k/divisions_per_bin
                CL_dw = 1.0*k/divisions_per_bin
                flag_1 = 1



        for k in range(central_value, upper_limit):
            value_int = real_Likelihood_Histogram.Integral(central_value - lower_limit, k - lower_limit + 1)
            if value_int > 0.475 and flag_2 == 0:    
                print " 95% confidence lower limit: ",1.0*k/divisions_per_bin
                CL_up = 1.0*k/divisions_per_bin
                flag_2 = 1

        for j in range(0, 100):
            marker_cl_2sigma.DrawMarker(CL_dw , 0.0001*j)
            marker_cl_2sigma.DrawMarker(CL_up , 0.0001*j)


        CL_up = (CL_up - CL_dw)/2
        print "Confidence Interval: ",CL_up
        print " Value of S: %f, \pm %f"%(position_x, CL_up)
        c1.SaveAs("normalized_Likelihood"+name+".pdf")


    return position_x
    


#______________________________________________________________

def ChiSquare_OneVariable_h2(lower_limit = 0, upper_limit = 10, divisions_per_bin = 1.0, doit = 0):

    ChiSquare_Histogram = TH1D("chi", "chi", upper_limit - lower_limit, lower_limit/divisions_per_bin, upper_limit/divisions_per_bin)

    chi_square = 999999.
    position_x = 0

    for x1 in range (lower_limit, upper_limit+1):
        xx = 1.0*x1/divisions_per_bin
        chi_square_temp = 0.
                
        for ii in range (1, 11):
            chi_square_temp = chi_square_temp + (1.0*xx - 1.0*h2.GetBinContent(ii))*(1.0*xx - 1.0*h2.GetBinContent(ii))/(1.0*xx)

        ChiSquare_Histogram.SetBinContent(x1-lower_limit, chi_square_temp)
        
        if  chi_square_temp < abs(chi_square):
            chi_square = chi_square_temp
            position_x = xx

    bin_err_up = 0.0
    bin_err_up_2sigma = 0.0
    bin_err_dw_2sigma = 0.0
    bin_err_dw = 0.0
    has_it_up = 0
    has_it_dw = 0
    has_it_up_2sigma = 0
    has_it_dw_2sigma = 0
    chi_square_up = 0.0
    chi_square_up_2sigma = 0.0
    chi_square_dw = 0.0
    chi_square_dw_2sigma = 0.0


    pre = "Long"


    if doit==1:
        pre = ""

        for x1 in range (position_x*divisions_per_bin, upper_limit+1):
            xx = 1.0*x1/divisions_per_bin
            chi_square_temp = 0.
            for ii in range (1, 11):
                chi_square_temp = chi_square_temp + (1.0*xx - 1.0*h2.GetBinContent(ii))*(1.0*xx - 1.0*h2.GetBinContent(ii))/(1.0*xx)

            if abs(chi_square_temp)> abs(chi_square)+thresholdchisquared and has_it_up == 0:
                has_it_up = 1
                bin_err_up = xx
                chi_square_up = chi_square_temp

            if abs(chi_square_temp)> abs(chi_square)+thresholdchisquared+1  and has_it_up_2sigma == 0:
                has_it_up_2sigma = 1
                bin_err_up_2sigma = xx
                chi_square_up_2sigma = chi_square_temp     
               
            
        for x1 in range (lower_limit, position_x*divisions_per_bin):
            xx = 1.0*x1/divisions_per_bin
            chi_square_temp = 0.
            for ii in range (1, 11):
                chi_square_temp = chi_square_temp + (1.0*xx - 1.0*h2.GetBinContent(ii))*(1.0*xx - 1.0*h2.GetBinContent(ii))/(1.0*xx)


            if abs(chi_square_temp)< abs(chi_square)+thresholdchisquared and has_it_dw == 0:
                has_it_dw = 1
                bin_err_dw = xx
                chi_square_dw = chi_square_temp
                
            if abs(chi_square_temp)< abs(chi_square)+thresholdchisquared + 1 and has_it_dw_2sigma == 0:
                has_it_dw_2sigma = 1
                bin_err_dw_2sigma = xx
                chi_square_dw_2sigma = chi_square_temp
                
        print "Chi Square: ",chi_square
        print "X position: ",position_x
        print "Error up x: ",bin_err_up
        print "Chi Square err up: ",chi_square_up
        print "Error dw x: ",bin_err_dw
        print "Chi Square err dw: ",chi_square_dw



        ChiSquare_Histogram.SetTitle("Calculating A from \Chi^{2} for h2")
        ChiSquare_Histogram.GetYaxis().SetTitle(" \Chi^{2} ")
        ChiSquare_Histogram.GetXaxis().SetTitle(" A [Data] ")
        ChiSquare_Histogram.SetStats(0)
        ChiSquare_Histogram.SetLineColor(kOrange)
        ChiSquare_Histogram.SetLineWidth(4)
        ChiSquare_Histogram.Draw("")

        
        text2 = TLatex(position_x, 2*chi_square, " A = %.2f "%position_x)
        marker_position.DrawMarker(position_x, chi_square)
        text2.Draw("same")


        marker_1sigma.DrawMarker(bin_err_up, chi_square_up)
        marker_1sigma.DrawMarker(bin_err_dw, chi_square_dw)


        marker_2sigma.DrawMarker(bin_err_up_2sigma, chi_square_up_2sigma)
        marker_2sigma.DrawMarker(bin_err_dw_2sigma, chi_square_dw_2sigma)



        c1.SaveAs("Chi_Square_h2.pdf")

    return position_x
#__________________________________________________________________________________________________
def ChiSquare_OneVariable_h3(lower_limit = 0, upper_limit = 10, divisions_per_bin = 1.0, doit = 0):

    ChiSquare_Histogram = TH1D("chi", "chi", upper_limit - lower_limit, lower_limit/divisions_per_bin, upper_limit/divisions_per_bin)

    chi_square = 999999.
    position_x = 0

    for x1 in range (lower_limit, upper_limit+1):
        xx = 1.0*x1/divisions_per_bin
        chi_square_temp = 0.
                
        for ii in range (1, 11):
            chi_square_temp = chi_square_temp + (1.0*xx - 1.0*h3.GetBinContent(ii))*(1.0*xx - 1.0*h3.GetBinContent(ii))/(1.0*xx)

        ChiSquare_Histogram.SetBinContent(x1-lower_limit, chi_square_temp)
        
        if  chi_square_temp < abs(chi_square):
            chi_square = chi_square_temp
            position_x = xx


    bin_err_up = 0.0
    bin_err_up_2sigma = 0.0
    bin_err_dw_2sigma = 0.0
    bin_err_dw = 0.0
    has_it_up = 0
    has_it_dw = 0
    has_it_up_2sigma = 0
    has_it_dw_2sigma = 0
    chi_square_up = 0.0
    chi_square_up_2sigma = 0.0
    chi_square_dw = 0.0
    chi_square_dw_2sigma = 0.0


    pre = "Long"


    if doit==1:
        pre = ""

        for x1 in range (position_x*divisions_per_bin, upper_limit+1):
            xx = 1.0*x1/divisions_per_bin
            chi_square_temp = 0.
            for ii in range (1, 11):
                chi_square_temp = chi_square_temp + (1.0*xx - 1.0*h3.GetBinContent(ii))*(1.0*xx - 1.0*h3.GetBinContent(ii))/(1.0*xx)

            if abs(chi_square_temp)> abs(chi_square)+thresholdchisquared and has_it_up == 0:
                has_it_up = 1
                bin_err_up = xx
                chi_square_up = chi_square_temp

            if abs(chi_square_temp)> abs(chi_square)+thresholdchisquared+1  and has_it_up_2sigma == 0:
                has_it_up_2sigma = 1
                bin_err_up_2sigma = xx
                chi_square_up_2sigma = chi_square_temp     
               
            
        for x1 in range (lower_limit, position_x*divisions_per_bin):
            xx = 1.0*x1/divisions_per_bin
            chi_square_temp = 0.
            for ii in range (1, 11):
                chi_square_temp = chi_square_temp + (1.0*xx - 1.0*h3.GetBinContent(ii))*(1.0*xx - 1.0*h3.GetBinContent(ii))/(1.0*xx)


            if abs(chi_square_temp)< abs(chi_square)+thresholdchisquared and has_it_dw == 0:
                has_it_dw = 1
                bin_err_dw = xx
                chi_square_dw = chi_square_temp
                
            if abs(chi_square_temp)< abs(chi_square)+thresholdchisquared + 1 and has_it_dw_2sigma == 0:
                has_it_dw_2sigma = 1
                bin_err_dw_2sigma = xx
                chi_square_dw_2sigma = chi_square_temp
                
    print "Chi Square: ",chi_square
    print "X position: ",position_x
    print "Error up x: ",bin_err_up
    print "Chi Square err up: ",chi_square_up
    print "Error dw x: ",bin_err_dw
    print "Chi Square err dw: ",chi_square_dw


    
       
    ChiSquare_Histogram.SetTitle("Calculating A from \Chi^{2} for h3")
    ChiSquare_Histogram.GetYaxis().SetTitle(" \Chi^{2} ")
    ChiSquare_Histogram.GetXaxis().SetTitle(" A [Data] ")
    ChiSquare_Histogram.SetStats(0)
    ChiSquare_Histogram.SetLineColor(kOrange)
    ChiSquare_Histogram.SetLineWidth(4)
    ChiSquare_Histogram.Draw("")

    text2 = TLatex(0.96*position_x, 3*chi_square, " A = %.2f "%position_x)
    marker_position.DrawMarker(position_x, chi_square)
    text2.Draw("same")


    marker_1sigma.DrawMarker(bin_err_up, chi_square_up)
    marker_1sigma.DrawMarker(bin_err_dw, chi_square_dw)


    marker_2sigma.DrawMarker(bin_err_up_2sigma, chi_square_up_2sigma)
    marker_2sigma.DrawMarker(bin_err_dw_2sigma, chi_square_dw_2sigma)



    c1.SaveAs("Chi_Square_h3.pdf")



#____________________________________________-
def Likelihood_OneVariable_h2(hx, lower_limit= 0, upper_limit = 10, divisions_per_bin= 1.0, doit=0):

    Likelihood_Histogram =TH1D( "Likelihood_Histogram", "Likelihood_Histogram", (upper_limit - lower_limit), lower_limit/divisions_per_bin, upper_limit/divisions_per_bin)#, 500, -200, 0)
    real_Likelihood_Histogram =TH1D( "real_Likelihood_Histogram", "real_Likelihood_Histogram", (upper_limit - lower_limit), lower_limit/divisions_per_bin, upper_limit/divisions_per_bin)#, 500, -200, 0)
    likelihood = - 99999.
    position_x = 0
    
    for x1 in range (lower_limit, upper_limit+1):
        xx = 1.0*x1/divisions_per_bin
        likelihood_temp = 0.0
        for ii in range (1, 11):
                f_value = xx
                likelihood_temp = likelihood_temp + (hx.GetBinContent(ii))*log(f_value) - f_value - log( Factorial( hx.GetBinContent(ii)))

        Likelihood_Histogram.SetBinContent(x1 - lower_limit, likelihood_temp)
        real_Likelihood_Histogram.SetBinContent(x1 - lower_limit, 1.0*exp(likelihood_temp))
        if abs(likelihood_temp) < abs(likelihood):    
            likelihood = likelihood_temp
            position_x = xx

    bin_err_up = 0.0
    bin_err_up_2sigma = 0.0
    bin_err_dw_2sigma = 0.0
    bin_err_dw = 0.0
    has_it_up = 0
    has_it_dw = 0
    has_it_up_2sigma = 0
    has_it_dw_2sigma = 0
    likelihood_up = 0.0
    likelihood_up_2sigma = 0.0
    likelihood_dw = 0.0
    likelihood_dw_2sigma = 0.0
    pre = "Long"

    if doit==1:
        pre = ""
        for x1 in range (position_x*divisions_per_bin, upper_limit):
            xx = 1.0*x1/divisions_per_bin
            likelihood_temp = 0.0
            for ii in range (1, 11):
                    f_value = xx
                    likelihood_temp = likelihood_temp + (hx.GetBinContent(ii))*log(f_value) - f_value - log( Factorial( hx.GetBinContent(ii)))
                   
            if abs(likelihood_temp)> abs(likelihood)+thresholdlikelihood and has_it_up == 0:
                has_it_up = 1
                bin_err_up = xx
                likelihood_up = likelihood_temp
            if abs(likelihood_temp)> abs(likelihood)+thresholdlikelihood +1  and has_it_up_2sigma == 0:
                has_it_up_2sigma = 1
                bin_err_up_2sigma = xx
                likelihood_up_2sigma = likelihood_temp     
        for x1 in range (lower_limit, position_x*divisions_per_bin):
            xx = 1.0*x1/divisions_per_bin
            likelihood_temp = 0.0
            for ii in range (1, 11):
                    f_value = xx
                    likelihood_temp = likelihood_temp + (hx.GetBinContent(ii))*log(f_value) - f_value - log( Factorial( hx.GetBinContent(ii)))

            if abs(likelihood_temp)< abs(likelihood)+thresholdlikelihood and has_it_dw == 0:
                has_it_dw = 1
                bin_err_dw = xx
                likelihood_dw = likelihood_temp
            if abs(likelihood_temp)< abs(likelihood)+thresholdlikelihood + 1and has_it_dw_2sigma == 0:
                has_it_dw_2sigma = 1
                bin_err_dw_2sigma = xx
                likelihood_dw_2sigma = likelihood_temp
                
    print "Likelihood: ",likelihood
    print "X position: ",position_x
    print "Error up x: ",bin_err_up
    print "Likelihood err up: ",likelihood_up
    print "Error dw x: ",bin_err_dw
    print "Likelihood err dw: ",likelihood_dw

    Likelihood_Histogram.SetTitle("Calculating A (h2) from Likelihold")
    Likelihood_Histogram.GetYaxis().SetTitle("log [ Likelihold]")
    Likelihood_Histogram.GetXaxis().SetTitle(" A  (h2) [Data] ")
    Likelihood_Histogram.SetStats(0)
    Likelihood_Histogram.SetLineColor(kOrange)
    Likelihood_Histogram.SetLineWidth(4)
    Likelihood_Histogram.Draw("")
    text = TLatex(0.95*position_x, 5*likelihood, " A = %.2f "%position_x)
    
    marker_position.DrawMarker(position_x, likelihood)
    text.Draw("same")

    marker_1sigma.DrawMarker(bin_err_up, likelihood_up)
    marker_1sigma.DrawMarker(bin_err_dw, likelihood_dw)
    c1.SaveAs(pre+"Likelihood_1var_h2.pdf")

    real_Likelihood_Histogram.SetTitle("Calculating A (h2) from Likelihold")
    real_Likelihood_Histogram.GetYaxis().SetTitle("Likelihold")
    real_Likelihood_Histogram.GetXaxis().SetTitle(" A (h2) [Data] ")
    real_Likelihood_Histogram.SetStats(0)
    real_Likelihood_Histogram.SetLineColor(kOrange)
    real_Likelihood_Histogram.SetLineWidth(4)
    real_Likelihood_Histogram.Draw("")
    text2 = TLatex(position_x, exp(likelihood), " A = %.2f "%position_x)
    marker_position.DrawMarker(position_x, exp(likelihood))
    marker_1sigma.DrawMarker(bin_err_up, exp(likelihood_up))
    marker_1sigma.DrawMarker(bin_err_dw, exp(likelihood_dw))
    marker_2sigma.DrawMarker(bin_err_up_2sigma, exp(likelihood_up_2sigma))
    marker_2sigma.DrawMarker(bin_err_dw_2sigma, exp(likelihood_dw_2sigma))
 
    text2.Draw("same")

    c1.SaveAs(pre+"real_Likelihood_1var_h2.pdf")

    total_integral = real_Likelihood_Histogram.Integral(0, 999999999)
    real_Likelihood_Histogram.Scale(1/total_integral)
    
    real_Likelihood_Histogram.SetTitle("Calculating A (h2) from Normalized Likelihood")
    real_Likelihood_Histogram.GetYaxis().SetTitle("Likelihold (Normalized) ")
    real_Likelihood_Histogram.GetXaxis().SetTitle(" A (h2) [Data] ")
    axis = TGaxis()
    axis.SetMaxDigits(2)
    
    real_Likelihood_Histogram.SetStats(0)
    real_Likelihood_Histogram.SetLineColor(kOrange)
    real_Likelihood_Histogram.SetLineWidth(4)
    real_Likelihood_Histogram.Draw("")
    text2 = TLatex(position_x, exp(likelihood)*(1/total_integral), " A = %.2f "%position_x)


    marker_position.DrawMarker(position_x, exp(likelihood)*(1/total_integral))
    marker_1sigma.DrawMarker(bin_err_up, exp(likelihood_up)*(1/total_integral))
    marker_1sigma.DrawMarker(bin_err_dw, exp(likelihood_dw)*(1/total_integral))
    marker_2sigma.DrawMarker(bin_err_up_2sigma, exp(likelihood_up_2sigma)*(1/total_integral))
    marker_2sigma.DrawMarker(bin_err_dw_2sigma, exp(likelihood_dw_2sigma)*(1/total_integral))
    text2.Draw("same")


    total_integral = real_Likelihood_Histogram.Integral(0, 9999999999)
    print "Total: ",total_integral
    central_value = int(position_x*divisions_per_bin)
    flag_1 = 0
    flag_2 = 0
    CL_dw = 0.0
    CL_up = 0.0
    
    for k in range(lower_limit, central_value):
        value_int = real_Likelihood_Histogram.Integral(1, k - lower_limit + 1)
        if value_int > 0.025 and flag_1 == 0:    
            print " 95% confidence lower limit: ",1.0*k/divisions_per_bin
            CL_dw = 1.0*k/divisions_per_bin
            flag_1 = 1



    for k in range(central_value, upper_limit):
        value_int = real_Likelihood_Histogram.Integral(central_value - lower_limit, k - lower_limit + 1)
        if value_int > 0.475 and flag_2 == 0:    
            print " 95% confidence lower limit: ",1.0*k/divisions_per_bin
            CL_up = 1.0*k/divisions_per_bin
            flag_2 = 1

    for j in range(0, 400):
        marker_cl_2sigma.DrawMarker(CL_dw , 0.00001*j)
        marker_cl_2sigma.DrawMarker(CL_up , 0.00001*j)


    CL_up = (CL_up - CL_dw)/2
    print "Confidence Interval: ",CL_up
    print " Value of A: %f, \pm %f"%(position_x, CL_up)
    c1.SaveAs("normalized_Likelihood_h2.pdf")



    
    
    
################### PART 1 ############################33




#ChiSquare_TwoVariables(380, 500, 10, 250, 400, 10, 1)
Likelihood_TwoVariables(380, 500, 10, 250, 400, 10, 1)
#Likelihood_OneVariable(200, 500, 10, 1)



################### PART 2 ##############################33


#ChiSquare_OneVariable_h3(700, 1500, 100, 1)
#ChiSquare_OneVariable_h2(600, 2000, 100, 1)
#Likelihood_OneVariable(h2,0, 1600, 100, 1)


# Use is constant = 1 for h2 and h3. Use is_constant = 0 for h
#Pseudo_Experiment(1000, h, 0)
