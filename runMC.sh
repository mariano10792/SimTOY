#!/bin/bash
clear
echo "Entro al directorio /home/dario/CCD/toy "
cd /home/mariano/Desktop/SimToy
echo "Paso a ROOT 6 "
. /home/mariano/root6/bin/thisroot.sh
echo "Compilo toyMC "
make &&
echo "Corro loop_MC "
echo ""
python loop_MC_grilla.py 
echo ""
echo "Termine el loop de Monte Carlo " 
echo "Paso a ROOT 5" 
. /home/mariano/root5/bin/thisroot.sh 
echo "Entro al directorio /MC" 
cd MC/ 
echo "Corro loop_clustering " 
echo ""
python loop_clustering.py 
echo ""
echo "Sumo usando hadd los diferentes RUNs hechos con toyMC" 
python hadd.py 
#echo "Elimino los fits, pdf y root de más."
#python rm_root.py
#echo "Chau roots."
#python rm_pdf.py
#echo "Chau pdf."
#python rm_fits.py
#echo "Chau fits."

#echo "Muevo los Out* a la carpeta compare"
#mv Out* /home/dario/CCD/toy/compare
#echo "Entro al directorio /home/dario/CCD/toy/compare" 
#cd /home/dario/CCD/toy/compare 
#echo "Paso a ROOT 6" 
#. /home/dario/root6/bin/thisroot.sh 
#echo "Ejecuto MakePlots para comparar el Montecarlo con el experimento"
#root -l MakePlots.C
