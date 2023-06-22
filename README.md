#phase2_efficiencyStudies

mkdir < workArea >  
cd < workArea >  
cmsrel CMSSW_10_6_27   
cd CMSSW_10_6_27/src   
cmsenv  
git clone git@github.com:MaryHadley/phase2_efficiencyStudies.git  
cd phase2_efficiencyStudies  
mv * ..   
cd ..   
rm -rf phase2_efficiencyStudies  

#To run code   
root -l -b -q 'loop_ZplusY_effStudies.C++("< fileName.root >")'

root -l -b -q 'debug_loop_ZplusY_effStudies.C++("< fileName.root >")'  

