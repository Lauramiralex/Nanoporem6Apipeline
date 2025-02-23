Snakefile.py:
Hardgecoded wie Experimente heißen, benutzt 5.1.0
Produziert pileup.bed

Snakefile2.py:
Wildcards für experiment Namen, aber Inhalt der Wildcards oben unter experiments hardgecoded festgelegt in Array
Produziert pileup.bed

Snakefile3.py
Wildcards für experiment Namen, aber Inhalt der Wildcards oben unter experiments hardgecoded festgelegt in Array
Produziert pileup.bed
Produziert logs für jede einzelne Rule

Snakefile4.py
Wildcards für experiment, configfile festgelegt über config/config.yaml. config.yaml verweist auf sample_sheet.tsv Tabelle in welcher man Namen der Experimente hinterlegen muss zur Verarbeitung.
Wenn man die Control auch verwertet benötigt, muss man control files ebenfalls in experiment spalte eintragen
Produziert pileup.bed
Produziert logs für jede einzelne Rule

How to start:
Als erstes navigieren zu Ordner in dem SNakefile liegt.
Über Console mit "conda activate laurasenv" Environment öffnen.

Für Snakefile1-3 benötigt man Basecall Modell, Referenz Genom, pod5 Ordner und Snakefile im selben Ordner:
snakemake --cores all -s snakefile{}.py --jobs 1 

Für Snakefile4 benötigt man Basecall Modell, Referenz Genom, pod5 Ordner, Config Ordner und Snakefile im selben Ordner:
snakemake --cores all -s snakefile{}.py --jobs 1 
