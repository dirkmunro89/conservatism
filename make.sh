rm itr_*.eps
rm tmp1_*.tex
rm tmp2_*.tex

python svanb.py
pdflatex iters_svanb.tex

rm itr_*.eps
rm tmp1_*.tex
rm tmp2_*.tex

python groen.py
pdflatex iters_groen.tex
