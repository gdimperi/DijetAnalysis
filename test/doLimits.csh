#! /bin/csh
combine -M Asymptotic RS1000_datacard.txt --rMax=5   -m 1000 -n Normal
combine -M Asymptotic RS1500_datacard.txt --rMax=20  -m 1500 -n Normal
combine -M Asymptotic RS2000_datacard.txt --rMax=50  -m 2000 -n Normal
combine -M Asymptotic RS2500_datacard.txt --rMax=150 -m 2500 -n Normal
combine -M Asymptotic RS3000_datacard.txt --rMax=300 -m 3000 -n Normal

combine -M Asymptotic RS1000_sub_datacard.txt --rMax=1   -m 1000 -n Sub
combine -M Asymptotic RS1500_sub_datacard.txt --rMax=5   -m 1500 -n Sub
combine -M Asymptotic RS2000_sub_datacard.txt --rMax=10  -m 2000 -n Sub
combine -M Asymptotic RS2500_sub_datacard.txt --rMax=50  -m 2500 -n Sub
combine -M Asymptotic RS3000_sub_datacard.txt --rMax=200 -m 3000 -n Sub
