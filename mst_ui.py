from Tkinter import * 
from mst_logic import *
import tkMessageBox
import os, platform

WINDOW_HEIGHT = 350
WINDOW_WIDTH = 900

global startTime
global timeA
global timeB
global temp_choice

global fileInput # input/input.txt

global folderName # folder name to store in output
global fileAll # text file contains all galaxy info (note that z is z in this one, not C*z)
global fileGrp # text file contains all galaxy in group info (note that z is C*z in this one, not z)
global fileGrpCent # text file contains all group center info (note that z is C*z in this one, not z)
global C # speed of light (km/s) - use for d = C*z/100 (Mpc/h)
global trimTree # trim length (Mpc/h) - in process_tree()
global cutTree # minimum number of member required to be considered filament - in process_tree()
global scoopR # radius of scooper from node of filament (Mpc/h) - in scooper()
global scoopD # distance of scooper from edge of filament (Mpc/h) - in scooper()
global trimTendril # trim length (Mpc/h) - in process_tendril()
global cutTendril # minimum number of member required to be considered tendril - in process_tendril()
D = math.acos(-1.0) / 180 # acos(-1) = Pi


class MyApp:
	def __init__(self, parent):
		operating_sys = platform.system()
		print operating_sys
		if operating_sys == 'Windows':
			self.build_ui(parent,labelPad=0)
		else:
			self.build_ui(parent,labelPad=4)


	def build_ui(self, parent, labelPad=0):
		t_label1 = "output folder name: "
		t_label2 = "galaxy list file name (needs to be in the input folder, ex: gal.txt): "
		t_label3 = "galaxy in group list file name (needs to be in the input folder, ex: grpMEM.txt: "
		t_label4 = "group center list file name (needs to be in the input folder, ex: grpCENT.txt): "
		t_label5 = "speed of light (km/s): "
		t_label6 = "trim length - in process_tree() (Mpc/h): "
		t_label7 = "minimum number of member required to be considered filament: "
		t_label8 = "radius of scooper from node of filament (Mpc/h): "
		t_label9 = "distance of scooper from edge of filament (Mpc/h): "
		t_label10 = "trim length - in process_tendril(): "
		t_label11 = "minimum number of member required to be considered tendril: "

		labelFrame = Frame(parent)
		buttonFrame = Frame(parent)
		entryFrame = Frame(parent, bg='blue')


		Label(labelFrame, text=t_label1, pady=labelPad, padx=labelPad).pack(anchor = E)
		Label(labelFrame, text=t_label2, pady=labelPad, padx=labelPad).pack(anchor = E)
		Label(labelFrame, text=t_label3, pady=labelPad, padx=labelPad).pack(anchor = E)
		Label(labelFrame, text=t_label4, pady=labelPad, padx=labelPad).pack(anchor = E)
		Label(labelFrame, text=t_label5, pady=labelPad, padx=labelPad).pack(anchor = E)
		Label(labelFrame, text=t_label6, pady=labelPad, padx=labelPad).pack(anchor = E)
		Label(labelFrame, text=t_label7, pady=labelPad, padx=labelPad).pack(anchor = E)
		Label(labelFrame, text=t_label8, pady=labelPad, padx=labelPad).pack(anchor = E)
		Label(labelFrame, text=t_label9, pady=labelPad, padx=labelPad).pack(anchor = E)
		Label(labelFrame, text=t_label10, pady=labelPad, padx=labelPad).pack(anchor = E)
		Label(labelFrame, text=t_label11, pady=labelPad, padx=labelPad).pack(anchor = E)

		self.e_1=StringVar()
		self.e_2=StringVar()
		self.e_3=StringVar()
		self.e_4=StringVar()
		self.e_5=StringVar()
		self.e_6=StringVar()
		self.e_7=StringVar()
		self.e_8=StringVar()
		self.e_9=StringVar()
		self.e_10=StringVar()
		self.e_11=StringVar()

		self.e_1.set('LATESTchangeADDimprovedINPUTZSDSS(trimmedToOneSide,z<0.06,M<-18.5,halomass1>12,no.grpmem>=2)')
		self.e_2.set('gal.txt')
		self.e_3.set('grpMEM.txt')
		self.e_4.set('grpCENT.txt')
		self.e_5.set('299792.458')
		self.e_6.set('4')
		self.e_7.set('3')
		self.e_8.set('5')
		self.e_9.set('5')
		self.e_10.set('1')
		self.e_11.set('2')

		self.entry1 = Entry(entryFrame, textvariable=self.e_1).pack()
		self.entry2 = Entry(entryFrame, textvariable=self.e_2).pack()
		self.entry3 = Entry(entryFrame, textvariable=self.e_3).pack()
		self.entry4 = Entry(entryFrame, textvariable=self.e_4).pack()
		self.entry5 = Entry(entryFrame, textvariable=self.e_5, text='299792.458').pack()
		self.entry6 = Entry(entryFrame, textvariable=self.e_6).pack()
		self.entry7 = Entry(entryFrame, textvariable=self.e_7).pack()
		self.entry8 = Entry(entryFrame, textvariable=self.e_8).pack()
		self.entry9 = Entry(entryFrame, textvariable=self.e_9).pack()
		self.entry10 = Entry(entryFrame, textvariable=self.e_10).pack()
		self.entry11 = Entry(entryFrame, textvariable=self.e_11).pack()

		Button(buttonFrame, text='Build input file', command=self.click_build_input).pack()
		Button(buttonFrame, text='Run', command=self.click_run).pack()
		Label(buttonFrame, text='Check the console for output').pack()

		labelFrame.pack(side='left')
		entryFrame.pack(side='left')
		buttonFrame.pack(side='left')


	def click_build_input(self):
		if(not self.is_input_ready()):
			print 'input not ready'
		else:
			data = self.e_1.get() + '\n'
			data += self.e_2.get() + '\n'
			data += self.e_3.get() + '\n'
			data += self.e_4.get() + '\n'
			data += self.e_5.get() + '\n'
			data += self.e_6.get() + '\n'
			data += self.e_7.get() + '\n'
			data += self.e_8.get() + '\n'
			data += self.e_9.get() + '\n'
			data += self.e_10.get() + '\n'
			data += self.e_11.get()

			with open('input/input.txt', 'w+') as input_file:
				input_file.write(data)



	def is_input_ready(self):
		string_vars = [self.e_1,self.e_2,self.e_3,self.e_4,self.e_5,self.e_6,self.e_7,self.e_8,self.e_9,self.e_10,self.e_11]
		for s in string_vars:
			if not s.get() or s.get() == '':
				return False
		return True

	def click_run(self):
		if self.is_input_ready():
			self.run_logic()
		else:
			print 'input not ready'


	def run_logic(self):
		run_logic()


root = Tk()
root.title('MST')
root.geometry(str(WINDOW_WIDTH) + 'x' + str(WINDOW_HEIGHT))
myapp = MyApp(root)
root.mainloop()
