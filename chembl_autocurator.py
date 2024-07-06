import tkinter as tk
from tkinter import *
from tkinter import messagebox
from tkinter import filedialog
from tkinter import ttk
import os
from tkinter.filedialog import askopenfilename
import pandas as pd
from rdkit import Chem
from rdkit.Chem import SaltRemover

form = tk.Tk()
form.title("ChEMBL_autocurator")
form.geometry("650x380")

tab_parent = ttk.Notebook(form)

tab1 = tk.Frame(tab_parent, background="#00ffff")
tab_parent.add(tab1, text="Data preparation")

initialdir=os.getcwd()

ls=['Molecule ChEMBL ID','Molecule Name', 'Molecule Max Phase','Molecular Weight', '#RO5 Violations', 'AlogP','Smiles', 'Standard Type',
'Standard Relation','Standard Value', 'Standard Units','Assay ChEMBL ID', 'Assay Description', 'Assay Type', 'BAO Format ID', 'BAO Label', 'Assay Organism',
'Assay Cell Type', 'Target ChEMBL ID','Target Name','Target Organism','Target Type', 'Document ChEMBL ID']

def datatr():
    global filename1
    filename1 = askopenfilename(initialdir=initialdir,title = "Select downloaded ChEMBL csv file")
    firstEntryTabOne.delete(0, END)
    firstEntryTabOne.insert(0, filename1)
    global c_
    c_,d_=os.path.splitext(filename1)
    global file1
    dm=Criterion0.get()
    print(dm)
    if dm=="com":
       file1 = pd.read_csv(filename1)
    elif dm=="tab":
       file1 = pd.read_csv(filename1, delimiter='\t')
    elif dm=='sem':
       file1 = pd.read_csv(filename1, delimiter=';')


def smitosmi(smi):
    mol=Chem.MolFromSmiles(smi)
    remover = SaltRemover.SaltRemover();
    m = remover.StripMol(mol)
    canSmi= Chem.MolToSmiles(m)
    if len(canSmi.split('.'))>1:
       f=len(canSmi.split('.')[0])
       l=len(canSmi.split('.')[1])
       if f>l:
          canSmi=canSmi.split('.')[0]
       else:
          canSmi=canSmi.split('.')[1]
    elif len(canSmi.split('.'))==1:
         pass
    return canSmi


def file_process():
    df=file1
    df=df.fillna('NA')
    df=df[ls]
    df1_1=df[(df['Standard Relation']!='NA')]
    df1=df1_1[(df1_1['Standard Value']!='NA')]
    df2=df1[(df1['Standard Value']>0)]
    #df2.to_csv('df2.csv')
    df2=df2[df2['Smiles']!='NA']
    df2['Canonical SMILES']=df2.apply(lambda x:smitosmi(x['Smiles']),axis=1)
    x=df2['Standard Units'].unique().tolist()
    if 'ug.mL-1' in x:
       dfu=df2[df2['Standard Units']=='ug.mL-1']
       dfn=df2[df2['Standard Units']=='nM']
       dfu['Standard Value']=dfu.apply(lambda x: (x['Standard Value']/x['Molecular Weight'])*1000*1000, axis=1)
       dfu['Standard Units']='nM (converted from ug.mL-1)'
       df2=pd.concat([dfn,dfu], axis=0)
    else:
       df2_1=df2[df2['Standard Units']=='nM']
       df2=df2_1
       
       
    df2.to_csv('df2.csv', index=False)
    bao=Criterion4.get()
    if bao=='single protein format':
       df3=df2[df2['BAO Label']=='single protein format']
    elif bao=='cell-based format':
         df3=df2[df2['BAO Label']=='cell-based format']
    elif bao=='assay format':
         df3=df2[df2['BAO Label']=='assay format']
    elif bao=='all':
         df3=df2
    if df3.shape[0]==0:
       messagebox.showinfo('# Error: ',"No such BAO Label is present"+"\n")

    at=Criterion3.get()
    if at=='asB':
       df4=df3[df3['Assay Type']=='B']
    elif at=='asF':
       df4=df3[df3['Assay Type']=='F']
    elif at=='asA1':
       df4=df3[df3['Assay Type']=='A']
    elif at=='asA':
         df4=df3
    print('df4 shape is {}'.format(df4.shape))
    if df4.shape[0]==0:
       messagebox.showinfo('# Error: ',"No such assay type is present"+"\n")
    

    ro=Criterion5.get()
    if ro=='ro0':
       df5=df4[df4['#RO5 Violations']=='0']
    elif ro=='ro1':
         df5=df4[(df4['#RO5 Violations']=='0') | (df4['#RO5 Violations']=='1')]
    elif ro=='ro2':
         df5=df4[(df4['#RO5 Violations']=='0') | (df4['#RO5 Violations']=='1') | (df4['#RO5 Violations']=='2')]
    elif ro=='roN':
         df5=df4
    #df5.to_csv('df5.csv', index=False)
    cut=int(N1B1_t1.get())
    #print(df.shape,df1.shape,df2.shape,df3.shape, df4.shape,df5.shape)
    df6=df5[df5['Standard Relation']!='NA']
    sls=list(df6['Standard Relation'].unique())
    try:
       eq=df6[df6['Standard Relation']=="'='"]
       eq['Active']=eq.apply(lambda x: 1 if x['Standard Value']<cut else 0, axis=1)
    except:
       pass
    try:
       gt=df6[df6['Standard Relation']=="'>'"]
       gtc=gt[gt['Standard Value']>=cut]
       gtc['Active']=0
    except:
       pass
    
    try:
       lt=df6[df6['Standard Relation']=="'<'"]
       ltc=lt[lt['Standard Value']<=cut]
       ltc['Active']=1
       print(ltc.shape)
    except:
       pass
    ct2=pd.concat([eq,gtc,ltc],axis=0)
    ct2.to_csv('ct2.csv', index=False)
    ls3=list(set(ct2['Canonical SMILES']))
    print(len(ls3))
    ls4=list(ct2['Canonical SMILES'])
    print(len(ls4))
    x,y=[],[]
    for i in ls4:
        if i not in x:
           x.append(i)
        else:
           y.append(i)
    print(len(x),len(y))
    ct2['DUPLICATE']=ct2.apply(lambda x:'Duplicate' if x['Canonical SMILES'] in y else 'Not Duplicate', axis=1)
    #ct2.to_csv('ct2.csv', index=False)
    ct3=ct2[ct2['DUPLICATE']=='Duplicate']
    #ct4=ct3[['Canonical SMILES']]
    dd=ct3.set_index('Canonical SMILES')
    ls5,ls6=[],[]
    for i in dd.index:
        if len(dd.loc[i]['Active'].unique())>1:
           ls5.append(i)
        else:
           ls6.append(i)

    if len(ls5)>0:
        messagebox.showinfo('# WARNING: ',"Duplicate compounds were obtained with conflicting results. These compounds will be saved as Duplicate_Conflict.csv "+"\n")
        duplWConf=ct2[ct2['Canonical SMILES'].isin(ls5)]
        try:
            duplWConf.to_csv(str(c_)+"_Duplicate_Conflict.csv", index=False)
        except PermissionError:
               messagebox.showinfo('# Error: ',"Please close the file with the same name."+"\n")
    ct4=ct2[ct2['DUPLICATE']=='Not Duplicate']
    #ct4.to_csv('Not_DUPLICATE.csv', index=False)
    ls7=[]
    #print(ls6)
    for i in set(ls6):
        ls7.append(ct2[ct2['Canonical SMILES'].isin([i])].iloc[0:1,:])
        ndd=pd.concat(ls7, axis=0)
    print('ls7 len is {}'.format(len(ls7)))
    if len(ls7)>0:
       fd=pd.concat([ndd,ct4], axis=0)
    else:
       fd=ct4
    try:
       fd.to_csv(str(c_)+'_Processed.csv', index=False)
       #messagebox.showinfo('# Relax: ',"Final files are saved."+"\n")
    except PermissionError:
           messagebox.showinfo('# Error: ',"Please close the file with the same name."+"\n")


firstLabelTabOne = tk.Label(tab1, text="Select downloaded ChEMBL file",font=("Helvetica", 12))
firstLabelTabOne.place(x=30,y=70)
firstEntryTabOne = tk.Entry(tab1,text='',width=40)
firstEntryTabOne.place(x=300,y=73)
b5=tk.Button(tab1,text='Browse', command=datatr,font=("Helvetica", 10))
b5.place(x=560,y=75)


Criterion_Label0 = ttk.Label(tab1, text="File delimiter",font=("Helvetica", 12),anchor=W, justify=LEFT)
Criterion0 = StringVar()
Criterion0.set('com')
Criterion_xc = ttk.Radiobutton(tab1, text='comma', variable=Criterion0, value='com')
Criterion_xt = ttk.Radiobutton(tab1, text='tab', variable=Criterion0, value='tab')
Criterion_xs = ttk.Radiobutton(tab1, text='Semicolon', variable=Criterion0, value='sem')
Criterion_Label0.place(x=150,y=25)
Criterion_xc.place(x=250,y=25)
Criterion_xt.place(x=330,y=25)
Criterion_xs.place(x=380,y=25)

N1B1_t1 = Label(tab1, text= 'Activity cut-off (in nM)',font=("Helvetica", 12))
N1B1_t1.place(x=150,y=110)
N1B1_t1=Entry(tab1)
N1B1_t1.place(x=320,y=110)


Criterion_Label3 = ttk.Label(tab1, text="Assay type",font=("Helvetica", 12),anchor=W, justify=LEFT)
Criterion3 = StringVar()
Criterion3.set('asA')
Criterion_b = ttk.Radiobutton(tab1, text='B', variable=Criterion3, value='asB')
Criterion_f = ttk.Radiobutton(tab1, text='F', variable=Criterion3, value='asF')
Criterion_a1 = ttk.Radiobutton(tab1, text='A', variable=Criterion3, value='asA1')
Criterion_a = ttk.Radiobutton(tab1, text='All', variable=Criterion3, value='asA')
Criterion_Label3.place(x=50,y=150)
Criterion_b.place(x=150,y=150)
Criterion_f.place(x=200,y=150)
Criterion_a1.place(x=250,y=150)
Criterion_a.place(x=300,y=150)

Criterion_Label4 = ttk.Label(tab1, text="BAO Label",font=("Helvetica", 12),anchor=W, justify=LEFT)
Criterion4 = StringVar()
Criterion4.set('all')
Criterion_acc3 = ttk.Radiobutton(tab1, text='Single Protein Format', variable=Criterion4, value='single protein format')
Criterion_roc3 = ttk.Radiobutton(tab1, text='Cell-based Format', variable=Criterion4, value='cell-based format')
Criterion_roc4 = ttk.Radiobutton(tab1, text='Assay Format', variable=Criterion4, value='assay format')
Criterion_roc5 = ttk.Radiobutton(tab1, text='All', variable=Criterion4, value='all')
Criterion_Label4.place(x=50,y=200)
Criterion_acc3.place(x=150,y=200)
Criterion_roc3.place(x=300,y=200)
Criterion_roc4.place(x=430,y=200)
Criterion_roc5.place(x=530,y=200)

Criterion_Label5 = ttk.Label(tab1, text="Maximum Ro5 violation",font=("Helvetica", 12),anchor=W, justify=LEFT)
Criterion5 = StringVar()
Criterion5.set('roN')
Criterion_ro0 = ttk.Radiobutton(tab1, text='0', variable=Criterion5, value='ro0')
Criterion_ro1 = ttk.Radiobutton(tab1, text='1', variable=Criterion5, value='ro1')
Criterion_ro2 = ttk.Radiobutton(tab1, text='2', variable=Criterion5, value='ro2')
Criterion_roN = ttk.Radiobutton(tab1, text='None', variable=Criterion5, value='roN')
Criterion_Label5.place(x=50,y=250)
Criterion_ro0.place(x=250,y=250)
Criterion_ro1.place(x=300,y=250)
Criterion_ro2.place(x=350,y=250)
Criterion_roN.place(x=400,y=250)

b2=Button(tab1, text='Generate datasets', command=file_process,bg="green",font=("Helvetica", 10),anchor=W, justify=LEFT)
b2.place(x=310,y=300)

tab_parent.pack(expand=1, fill='both')

form.mainloop()
