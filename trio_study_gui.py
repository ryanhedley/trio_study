import sys, csv
import tkinter as tk
from tkinter import filedialog as fd 

# Requires GenotypeCalls, BeadPoolManifest, code2genotype in module folder from https://github.com/Illumina/BeadArrayFiles
from module import GenotypeCalls, BeadPoolManifest, code2genotype

#SET CONSTANTS
BUTTONLEFTPOS = 50
BUTTONRIGHTPOS = 500
FILELABELPOS = 220
FILELABELWIDTH = 60

class Alleles:
    def __init__(self, locus, proband, father, mother):
        self.locus = locus
        self.proband = proband
        self.father = father
        self.mother = mother

    def is_biparental(self):

        test = (
            (self.proband == 'AB') and 
            ((self.father == 'AA' and self.mother == 'BB') or
            (self.father == 'BB' and self.mother == 'AA'))
        )

        return test

    def is_maternal(self):

        test = (
            (
            (self.proband == 'AA') and
            ((self.father == 'BB' and self.mother == 'AB') or
            (self.father == 'BB' and self.mother == 'AA'))
            ) or
            (self.proband == 'BB') and
            ((self.father == 'AA' and self.mother == 'AB') or
            (self.father == 'AA' and self.mother == 'BB'))
        )
        
        return test

    def is_paternal(self):

        test = (
            (
            (self.proband == 'AA') and
            ((self.father == 'AA' and self.mother == 'BB') or
            (self.father == 'AB' and self.mother == 'BB'))
            ) or
            (self.proband == 'BB') and
            ((self.father == 'BB' and self.mother == 'AA') or
            (self.father == 'AB' and self.mother == 'AA'))
        )
        
        return test


def print_out(genotype_filename, genotype, chrom, base):
    # print(genotype.locus, genotype.proband, genotype.father, genotype.mother, end=' ')
    with open(genotype_filename, 'a', newline='') as csvfile:
        fieldnames = ['Chr', 'Position', 'Locus', 'Proband', 'Mother', 'Father', 'Inheritance']
        spamwriter = csv.DictWriter(csvfile, fieldnames=fieldnames)

        if genotype.is_biparental():
            spamwriter.writerow({'Chr': chrom, 'Position': base, 'Locus': genotype.locus, 'Proband': genotype.proband, 'Mother': genotype.mother, 'Father': genotype.father, 'Inheritance': 'Biparental'})
        elif genotype.is_paternal():
            spamwriter.writerow({'Chr': chrom, 'Position': base, 'Locus': genotype.locus, 'Proband': genotype.proband, 'Mother': genotype.mother, 'Father': genotype.father, 'Inheritance': 'Paternal'})
        elif genotype.is_maternal():
            spamwriter.writerow({'Chr': chrom, 'Position': base, 'Locus': genotype.locus, 'Proband': genotype.proband, 'Mother': genotype.mother, 'Father': genotype.father, 'Inheritance': 'Maternal'})
    # else:
    #     print('Uninformative')


def inherit(genotype, count, bi_inh, pat_inh, mat_inh):
    
    count += 1

    if genotype.is_biparental():
        bi_inh += 1
    elif genotype.is_paternal():
        pat_inh += 1
    elif genotype.is_maternal():
        mat_inh += 1

    return count, bi_inh, pat_inh, mat_inh


def report(report_filename, count, bi_inh, pat_inh, mat_inh, chrom, base):
    informative = bi_inh + pat_inh + mat_inh
    
    with open(report_filename, 'w') as report_file:
        #report_file.write("hello")

        
    
    
        report_file.write('***** REPORT *****\n\n')
        report_file.write('Total SNPs =\t' + str(count) + '\n')
        report_file.write('Informative SNPs =\t' + str(informative) + '\n')
        report_file.write('Biparental =\t' + str(bi_inh) + '(' + str(round((100*(bi_inh/informative)),1)) + '%)\n')
        report_file.write('Paternal =\t' + str(pat_inh) + '(' + str(round((100*(pat_inh/informative)),1)) + '%)\n')
        report_file.write('Maternal =\t' + str(mat_inh) + '(' + str(round((100*(mat_inh/informative)),1)) + '%)\n')
        report_file.write('\n******************\n')




def run_analysis():

    global chrom_choice, message

    chrom_of_interest = chrom_choice.get()


    # Intialise new CSV file for genotyping data with headers
    with open(genotype_filename, 'w', newline='') as csvfile:
        fieldnames = ['Chr', 'Position', 'Locus', 'Proband', 'Mother', 'Father', 'Inheritance']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

        writer.writeheader()

    message = 'Opening ' + manifest_filename
    # send_message(message)
    print(message)

    names = BeadPoolManifest( manifest_filename ).names
    chroms = BeadPoolManifest( manifest_filename).chroms
    bases = BeadPoolManifest( manifest_filename).map_infos

    message = 'Opening ' + proband_filename
    # send_message(message)
    print(message)
    proband_genotypes = GenotypeCalls( proband_filename ).get_genotypes()

    message = 'Opening ' + father_filename
    # send_message(message)
    print(message)
    father_genotypes = GenotypeCalls( father_filename ).get_genotypes()

    message = 'Opening ' + mother_filename
    # send_message(message)
    print(message)
    mother_genotypes = GenotypeCalls( mother_filename ).get_genotypes()


    count = 0; bi_inh = 0; pat_inh = 0; mat_inh = 0 #Set variables for biparental, paternal and maternal inheritance


    for (chrom, locus, base, proband_genotype, mother_genotype, father_genotype) in zip( chroms, names, bases, proband_genotypes, mother_genotypes, father_genotypes ):
        if chrom == chrom_of_interest:
            genotype = Alleles(locus, code2genotype[proband_genotype], code2genotype[father_genotype], code2genotype[mother_genotype])
            print_out(genotype_filename, genotype, chrom, base)
            count, bi_inh, pat_inh, mat_inh = inherit(genotype, count, bi_inh, pat_inh, mat_inh)


    
    message = "Genotypes Written: " + genotype_filename
    print(message)
    send_message(message)

    report(report_filename, count, bi_inh, pat_inh, mat_inh, chrom, base)

    message = "Report Written: " + report_filename
    print(message)
    send_message(message)
    



def get_proband_file():
    global proband_filename
    proband_filename = fd.askopenfilename(title= "Select Proband GTC File", filetypes = (("gtc files","*.gtc"),("all files","*.*"))) 
    tk.Label(root, text=proband_filename, width=FILELABELWIDTH, font=('Helvetica', 8), anchor='e', bg='light goldenrod').place(x=FILELABELPOS, y=60)

def get_mother_file():
    global mother_filename
    mother_filename = fd.askopenfilename(title= "Select Mother GTC File", filetypes = (("gtc files","*.gtc"),("all files","*.*")))
    tk.Label(root, text=mother_filename, width=FILELABELWIDTH, font=('Helvetica', 8), anchor='e', bg='LightPink1').place(x=FILELABELPOS, y=100)

def get_father_file():
    global father_filename
    father_filename = fd.askopenfilename(title= "Select Father GTC File", filetypes = (("gtc files","*.gtc"),("all files","*.*")))
    tk.Label(root, text=father_filename, width=FILELABELWIDTH, font=('Helvetica', 8), anchor='e', bg='LightSkyBlue1').place(x=FILELABELPOS, y=140)

def manifest_file():
    global manifest_filename
    manifest_filename = fd.askopenfilename(title= "Select Manifest File", filetypes = (("bpm files","*.bpm"),("all files","*.*")))
    tk.Label(root, text=manifest_filename, width=FILELABELWIDTH, font=('Helvetica', 8), anchor='e', bg='gray81').place(x=FILELABELPOS, y=180)

def report_file():
    global report_filename
    report_filename = fd.asksaveasfilename(initialfile="trio_report.txt", title= "Save Report File", defaultextension = [("txt files","*.txt")], filetypes = (("csv files","*.txt"),("all files","*.*")))
    tk.Label(root, text=report_filename, width=FILELABELWIDTH, font=('Helvetica', 8), anchor='e', bg='gray81').place(x=FILELABELPOS, y=340)

def genotype_file():
    global genotype_filename
    genotype_filename = fd.asksaveasfilename(initialfile="trio_genotypes.csv", title= "Save Genotype File", defaultextension = [("csv files","*.csv")], filetypes = (("csv files","*.csv"),("all files","*.*")))
    tk.Label(root, text=genotype_filename, width=FILELABELWIDTH, font=('Helvetica', 8), anchor='e', bg='gray81').place(x=FILELABELPOS, y=380)

def send_message(message):
    global msg_box
    msg_box.insert(tk.END, message + '\n')
    msg_box.place(x=BUTTONLEFTPOS, y=450)

def close_app():

    root.destroy()

#SET VARIABLES
chrom_no = ''
chrlist = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12',
            '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y']
message = ''

    
root = tk.Tk()
root.geometry("700x700")
#root.config(bg='gray81')


#TITLE BAR
tk.Label(root, text="Trio Study Calculator", font="Verdana 16 bold").place(x=BUTTONLEFTPOS, y=10)

#PROBAND
tk.Button(root, text='Select Proband GTC File', width=20, bg='light goldenrod',
    command=get_proband_file).place(x=BUTTONLEFTPOS, y=60)

#MOTHER
tk.Button(root, text='Select Mother GTC File', width=20, bg='LightPink1',
    command=get_mother_file).place(x=BUTTONLEFTPOS, y=100)

#FATHER
tk.Button(root, text='Select Father GTC File', width=20, bg='LightSkyBlue1',
    command=get_father_file).place(x=BUTTONLEFTPOS, y=140)

#MANIFEST
tk.Button(root, text='Select Manifest BPM File', width=20, bg='gray81',
    command=manifest_file).place(x=BUTTONLEFTPOS, y=180)

#CHROMOSOME_DROPDOWN
chrom_choice = tk.StringVar(root)
chrom_choice.set(chrlist[0])
chrom_no = tk.OptionMenu(root, chrom_choice, *chrlist)
chrom_no.config(width=4, bg='gray81', font=('Helvetica', 10))
chrom_no.grid(row=5, column=2)
chrom_no.place(x=BUTTONLEFTPOS+30, y=230)
tk.Label(root, text="Chromosome to analyse", width=20, font=('Helvetica', 10)).place(x=BUTTONLEFTPOS+130, y=235)

#REPORT FILE SAVE
tk.Button(root, text='Report Save Location', width=20, bg='gray81',
    command=report_file).place(x=BUTTONLEFTPOS, y=340)

#GENOTYLE FILE SAVE
tk.Button(root, text='Genotypes Save Location', width=20, bg='gray81',
    command=genotype_file).place(x=BUTTONLEFTPOS, y=380)

#MESSAGE BOX
msg_box = tk.Text(root, height=5, width=70)
msg_box.place(x=BUTTONLEFTPOS, y=450)
send_message('Please select GTC, BPM and save file locations \nand chromosome to analyse above\n')

#RUN
tk.Button(root, text='Run Analysis', width=20, bg='lime green',
    command=run_analysis).place(x=BUTTONLEFTPOS, y=600)

#QUIT
tk.Button(root, text='Quit', width=20, bg='firebrick1',
    command=close_app).place(x=BUTTONRIGHTPOS, y=600)


root.mainloop()


