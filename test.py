import statNS

#=====================================================================
#Example 0: Compute f-mode frequencies given an array of central density

#define central density sequence in the SI unit kg/m^3
RhocSI = (statNS.c_double*19)(6.5e17,7e17,7.2e17,8e17,8.2e17,8.4e17,8.6e17,8.8e17,9e17,9.2e17,9.5e17,9.8e17,1e18,1.05e18,1.08e18,1.1e18,1.12e18,1.14e18,1.2e18)

#initiate an EoS_t type variable to store your EoS info
myEoS = statNS.EoS_t()

#initiate a CompactStar_t type array for storing the output results,
#you may register a large enough array to store results when using longer input.
Results = (statNS.CompactStar_t * 18)()

#load EOS file "APR.txt" to myEoS, be careful the b"xxx" format is Python dedicated, strange isn't it?
statNS.loadEoS(statNS.pointer(myEoS), b"./EoS_lib/APR.txt")

#call f-mode multi-thread computation, it will return the results to Results[],
#the last two numbers represent: 18 the length of RhocSI[], 4 the total threads for computation.
statNS.fmode_mt(Results, statNS.pointer(myEoS), RhocSI, 18, 4)

#print results to console
print("\nExample 0 output:")
for pf in range(0,18):
    print("Rhoc={:.5e}, M={:.4f}, R={:.4f}, freq={:.4f}, dmpTime={:.4f}".format(
        Results[pf].Rho,
        Results[pf].M,
        Results[pf].r,
        Results[pf].freq,
        Results[pf].dampTime))


#=====================================================================
#Example 1: Compute f-mode frequencies given an array of mass

#define mass array sequence in the unit of solar mass
massArray = (statNS.c_double*16)(1.25, 1.325, 1.4, 1.44 ,1.46, 1.48, 1.5, 1.55, 1.62, 1.64, 1.67, 1.7123, 1.7456, 1.8257, 1.90123, 1.95)

#This rhocArray[] is used to store the central density, it needs a proper size not less than massArray[]
rhocArray = (statNS.c_double*16)()

#This function compute the central density that corresponds to the mass array one by one,
#The central density will be written into rhocArray[],
#the last two numbers are, 16 the length fo massArray[], 6 the desired threads to use
statNS.M2Rhoc_Arr_fm(rhocArray, statNS.pointer(myEoS), massArray, 16, 6)

#Call the f-mode multi-thread computation with rhocArray[], see also in Example 0.
statNS.fmode_mt(Results, statNS.pointer(myEoS), rhocArray, 16, 6)

#print results to console
print("\nExample 1 output:")
for pf in range(0,16):
    print("Rhoc={:.5e}, M={:.4f}, R={:.4f}, freq={:.4f}, dmpTime={:.4f}".format(
        Results[pf].Rho,
        Results[pf].M,
        Results[pf].r,
        Results[pf].freq,
        Results[pf].dampTime))


#=====================================================================
#Example 2: Compute TOV equation given an array of central density

#Similar to the fmode_mt() in Example 0, 18 is the length or RhocSI[]
#There will be in total 6 parallel worker threads.
statNS.solveTOV_mt(Results, statNS.pointer(myEoS), RhocSI, 18, 6)
#print results to console
print("\nExample 2 output:")
for pf in range(0,16):
    print("Rhoc={:.5e}, M={:.4f}, R={:.4f}, I={:.4f}, Lambda={:.4f}".format(
        Results[pf].Rho,
        Results[pf].M,
        Results[pf].r,
        Results[pf].I,
        Results[pf].Lambda))


#=====================================================================
#Example 3: Compute TOV equation given an array of mass

#Similar usage as M2Rhoc_Arr_fm() in Example 1,
#NOTICE: 
#integrator for fmode-xx related computation is different from that for solveTOV-xx
#So I implement the central density determining function M2Rhoc() in two.
statNS.M2Rhoc_Arr_s(rhocArray, statNS.pointer(myEoS), massArray, 16, 4)
	
#Call the multi-thread TOV solver as Example 2
statNS.solveTOV_mt(Results, statNS.pointer(myEoS), rhocArray, 16, 4)
	
#print results to console
print("\nExample 3 output:")
for pf in range(0,16):
    print("Rhoc={:.5e}, M={:.4f}, R={:.4f}, I={:.4f}, Lambda={:.4f}".format(
        Results[pf].Rho,
        Results[pf].M,
        Results[pf].r,
        Results[pf].I,
        Results[pf].Lambda))

#=====================================================================
#Example 4: Compute TOV equation given an array of central density and EoS parameters

#We integrate the amEoS (Nai-Bo Zhang and Bao-An Li)
#As an example here we choose XL=58.7 MeV, Ksym=-105 MeV, Jsym=400 MeV, J0=115 MeV
if(statNS.genAmEoS(statNS.pointer(myEoS), 58.7, -105, 400, 115, None)<0):
    print("\nAmEoS failed!\n")
    exit
	
statNS.solveTOV_mt(Results, statNS.pointer(myEoS), RhocSI, 16, 4)
#print results to console
print("\nExample 4 output:")
for pf in range(0,16):
    print("Rhoc={:.5e}, M={:.4f}, R={:.4f}, I={:.4f}, Lambda={:.4f}".format(
        Results[pf].Rho,
        Results[pf].M,
        Results[pf].r,
        Results[pf].I,
        Results[pf].Lambda))


#=====================================================================
#Example 5: Use EoS_opt_t to alter the default density grid of amEoS

#Firstly, you need to clearify an EoS_opt_t object
myeosopt = statNS.EoS_opt_t()
statNS.set_EoS_default_opt(statNS.pointer(myeosopt))

#set the core segment's length to 250
#maximum hadronic density set to 8.0 saturation density
#disable the default dU (hadronic density's stepsize)
myeosopt.core_length=250
myeosopt.maxU=8.0
myeosopt.dU=-1

if(statNS.genAmEoS(statNS.pointer(myEoS), 58.7, -105, 400, 115, statNS.pointer(myeosopt))<0):
    print("\nAmEoS failed!\n")
    exit

statNS.solveTOV_mt(Results, statNS.pointer(myEoS), RhocSI, 16, 4)
#print results to console
print("\nExample 5 output:")
for pf in range(0,16):
    print("Rhoc={:.5e}, M={:.4f}, R={:.4f}, I={:.4f}, Lambda={:.4f}".format(
        Results[pf].Rho,
        Results[pf].M,
        Results[pf].r,
        Results[pf].I,
        Results[pf].Lambda))

#save the generated EoS to file
statNS.saveEoS(statNS.pointer(myEoS), b"EoS_lib/test_EoSgen.txt")
print("\nTest EoS saved to EoS_lib/test_EoSgen.txt\n")