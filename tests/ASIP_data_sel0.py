T_ASIP=[-3893,-2516,-638,-448,-116,0]
M_ASIP=[10,22,20,20,36,38]
I_ASIP=[0,1,15,12,15,18]

print " INPUT DATA: "
print 'Times_ASIP ',T_ASIP
print 'M_ASIP ',M_ASIP
print 'I_ASIP ',I_ASIP

I_=I_ASIP
M_=M_ASIP
T_=T_ASIP

#Upper_bounds=(-2516,59,100000)
#Lower_bounds=(-5000,-59,1000)

#Upper_bounds=(-2516,59,50000)
#Lower_bounds=(-5000,-59,5000)

fixed_params_=[None,0,1830]
#fixed_params_=[0,0,0]

#Upper_bounds=(-2516,59,50000)
#Lower_bounds=(-10000,-10,20000)

#Upper_bounds=(-2516,500,1e5)
#Lower_bounds=(-4000-10e5,-20,1e5)

#Upper_bounds=(-2516,500,5000)
#Lower_bounds=(-4000-10e4,-20,200)

Upper_bounds_=(-2516,40,5000)
Lower_bounds_=(-1e4,-40,200)

#Upper_bounds_exhaustive=(-2516,40,5000)
#Lower_bounds_exhaustive=(-8000,-40,200)

#Upper_bounds_boot_=(-2516,200,120000)
#Lower_bounds_boot_=(-8000,-200,200)


#print Lower_bounds
#Upper_bounds=(-2516,59,50000)
#Lower_bounds=(-10000,-10,10000)

dominance_=0.

# epistatic efect between ASIP and MC1R.  
# A-E- : bay , wild type
# aaE- : black, not wild.
# --ee : chestnut, not wild.

#round(-115.72-331.92-190.6-1878.04-1376.36-0)
