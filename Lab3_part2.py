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



#_________________________________________________________________________________
def Likelihood_second(A_0, hx = TH1F("ac", "ac", 10, 0, 10), is_constant = 0):

    nmax = hx.GetSize()-1
    sum2 = 0.0      
    for i in range (1, nmax):


        ff = "(%f)*exp((-1.0)*(%f)*(x)) + (%f)*exp((-1.0)*(x - %f)*(x-%f)/(%f*%f*2))"%(42.8, alpha, 0, mass_0, mass_0, sigma, sigma)
        function = TF1("function", ff, 0, 41)
        f_value = 1.0*function.Eval(1.0*max_bin_value*i/bins_number, 0, 0)

      
        sum2 = sum2 + (hx.GetBinContent(i))*log(f_value) - f_value - log( Factorial( hx.GetBinContent(i)))

        

    return sum2






#_____________________________________________________________________________________________________
def Pseudo_Experiment(number_of_experiments = 10, hx = TH1D("fakeasdfas", "fakeasfdasa", 50, 0, 5), is_constant = 1):


   
    A_0 = 0
    position_x = Likelihood_second(A_0, hx, is_constant)

    
    
    fake_experiment = TH1D("fake", "fake", 10, 0, 10)
    Likelihood_Histogram = TH1D("ll", "ll", 200, int(position_x*1.8), int(position_x*0.07))




    for j in range (0, number_of_experiments):

        binn = 10
        if is_constant == 0:
            binn = 50
    
        for k in range(0, binn):
            

            if is_constant == 0:
                 xx = TRandom1().PoissonD(mass_0)
                 ff = "(42.80)*exp((-1.0)*(%f)*(x)) + (31.1)*exp((-1.0)*(x - %f)*(x-%f)/(%f*%f*2))"%(alpha, mass_0, mass_0, sigma, sigma)
                 function = TF1("func", ff, 0, 50)


                 fake_experiment.SetBinContent(k+1, function.Eval(xx, 0,0))


        



        if is_constant == 0:

                Likelihood_Histogram.Fill(Likelihood_second(31.8, fake_experiment, is_constant))            

    
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


    total_integral = Likelihood_Histogram.Integral(1, 99999)


    partial_integral = Likelihood_Histogram.Integral(1, bin_x)/total_integral



    text3 = TLatex(1.6*position_x, 0.5*maximum_y, "NULL p value: %.3f"%partial_integral)    

    text3.Draw("same")
    
    
    c1.SaveAs("Nullfake"+name+".pdf")

    
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

 

    return position_x
    


    

################### PART 2 ##############################33


#ChiSquare_OneVariable_h3(700, 1500, 100, 1)
#ChiSquare_OneVariable_h2(600, 2000, 100, 1)
#Likelihood_OneVariable(h2,0, 1600, 100, 1)


# Use is constant = 1 for h2 and h3. Use is_constant = 0 for h
Pseudo_Experiment(10000, h, 0)
