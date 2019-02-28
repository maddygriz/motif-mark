#!/usr/bin/env python

import cairo
import re
import random
import argparse

#*******************************************************************************
#Function	: get_arguments
#Description: get the arguments from the command line for generalization 
#			  of the program
#Parameters	: none
#Returned	: parse_args - the arguments needed for the program
#******************************************************************************* 

def get_arguments():
	parser = argparse.ArgumentParser (description = "creating a figure for a gene and location of motifs")
	parser.add_argument("-g", "--genes", help="fasta file of genes", required=True, type=str)
	parser.add_argument("-m", "--motifs", help="motifs needed to look for, max = 10", required=True, type=str)
	parser.add_argument("-o", "--output", help="location for outputing images", required=True, type=str)
	parser.add_argument("-f", "--file_type", help="type of file for image files, default .svg, suggestions: .png, .jpeg", required=False, type=str, default=".svg")
	parser.add_argument("-c", "--colors", help="color palletes for images, options: general, blues, reds, random", required=False, type=str, default="general", choices=['general','blues','reds',"random"])
	return parser.parse_args()
args = get_arguments()

#*******************************************************************************
#Function 	: exoc_loc
#Description: locate exons in FASTA file (capital nucleotides) 
#Parameters : string - fasta sequence
#Returned 	: exon - dictionary with start and stop values for each exon
#*******************************************************************************

def exon_loc(string):
	first = True
	counter = 0
	exon_num = 0
	exons = {}

	for char in string:
		counter += 1
		if ord(char) >= 65 and ord(char) <= 90:
			if first == True:
				exon_num += 1
				start = counter
			else:
				end = counter
			first = False
		else:
			if first == False:
				exons[exon_num] = (start, end)
			first = True

	return(exons)
	
#*******************************************************************************
#Function 	: cairo_init
#Description: initiate cairo canvas
#Parameters : names - name of saved svg
#Returned 	: surface - part of canvas needed for Cairo
#			  context - part of canvas needed for Cairo
#*******************************************************************************	
	
def cairo_init(name, file_type):
	width, height = 800, 500
	if file_type == ".svg":
		#create the coordinates to display your graphic, desginate output
		surface = cairo.SVGSurface(name,width, height)
	else:
		surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, width, height)
	#create the coordinates you will be drawing on (like a transparency) - you can create a transformation matrix
	context = cairo.Context(surface)
	#context.scale(width,height) #will set your drawing surface to a 0.0-1.0 scale
	context.scale(width,height)
	
	return surface, context

#input fasta file of genes
#input text file with motifs, one per line
genes_txt = args.genes
motifs_txt = args.motifs
output_folder = args.output

#check for file types
if ".fa" not in genes_txt:
	exit("gene file was not a fasta file with extension .fa or .fasta")
if ".txt" not in motifs_txt:
	exit("motifs file did not have extention .txt")

#initailize motif dictionary 
#key = given motif 
#value = regex search 

#IUPAC codes https://www.bioinformatics.org/sms/iupac.html
nuc_codes = {}
nuc_codes["A"] = "[Aa]"
nuc_codes["T"] = "[Tt]"
nuc_codes["G"] = "[Gg]"
nuc_codes["C"] = "[Cc]"
nuc_codes["R"] = "[AaGg]"
nuc_codes["Y"] = "[CcTt]"
nuc_codes["S"] = "[GgCc]"
nuc_codes["W"] = "[AaTt]"
nuc_codes["K"] = "[GgTt]"
nuc_codes["M"] = "[AaCc]"
nuc_codes["B"] = "[CcGgTt]"
nuc_codes["D"] = "[AaGgTt]"
nuc_codes["H"] = "[AaCcTt]"
nuc_codes["V"] = "[AaCcGg]"
nuc_codes["N"] = "[AaCcGgTt]"

motifs = {}

with open(motifs_txt, "r") as fh:
	count = 0
	for line in fh:
		count += 1
		line = line.strip()
		motif = ""
		if count == 11:
			exit("Too many motifs for the system to handle")
		for char in line:
			if char.capitalize() in nuc_codes:
				motif += nuc_codes[char.capitalize()]
			else:
				exit("One of the given motifs contains an unidentified IUPAC nucleotide code")
		motifs[line] = motif

		
#initialize dictionary with genes
#key = gene names
#value = gene sequence

gene_info = {}
gene_len = {}

read = ""
first = True
with open(genes_txt, "r") as gene:
    for line in gene:
        line = line.strip()
        if line[0] != ">":
            read += line
        else:
            if first != True:
                gene_info[gene_name] = read
                read = ""
            else:
                first = False
            gene = re.search(">(\S+)", line)
            gene_name = gene.group(1)
            startLoc = re.search(".+ (.+:[0-9]+-[0-9]+)", line)
            gene_len[gene_name] = startLoc.group(1)
                
gene_info[gene_name] = read

#set colors palletes for each motif
if args.colors == "general":
	colors = (.9,.6,0), (.8,0.1,0.1), (.6,0.2,0.8), (0.1,0.9,0.9), (.4, 0.9, 0.2), (.0, 0.0, 0.9), (.9, 0.9,0), (0.9,0.3,0.8), (.5,0.5,0.5), (.3,0.3,0.9)
if args.colors == "blues":
	colors = (.5,1,1), (0,.2,1), (.4,.8,.6), (0,.6,.8), (.4,.8,1), (0,1,.2), (.4,.9,.7), (0,.8,.5), (.4,.4,.8), (0,.2,1)
if args.colors == "reds":
	colors = (1,0,0), (1,0,.5), (1,.4,1), (.8,.2,.5), (.6,.2,.9), (1,.4,0), (.6,.4,.7), (0.9,.5,.5), (.7,.4,1), (1,.4,.6)
if args.colors == "random":
	colors = []
	for i in range(1,11):
		red = random.uniform(0,1)
		green = random.uniform(0,1)
		blue = random.uniform(0,1)
		color = red, green, blue
		colors.append(color)

for gene in gene_info:    
	#cairo drawing gene length
	surface, context = cairo_init(output_folder+gene+args.file_type, args.file_type)
	context.set_line_width(.005)
	context.move_to(.05,.7)        #(x,y)
	context.line_to(.95,.7)
	context.stroke()
	
	#scaling and shift for equal image size and location
	#multiplied by .9 to account for length of gene line
	#motif counter initialized - negative one for zero-based list
	scale = (1/len(gene_info[gene])) * .9
	shift = .05
	key_count = -1
	
	#add exon location to image
	exons = exon_loc(gene_info[gene])
	for exon in exons:
		context.rectangle((exons[exon][0]*scale) + shift,.65,(exons[exon][1]-exons[exon][0])*scale,.1)        #(x0,y0,width,height)
		context.fill()
	
	#add gene name to image
	context.select_font_face("Purisa", cairo.FONT_SLANT_NORMAL, 
		cairo.FONT_WEIGHT_NORMAL)
	context.set_font_size(.1)
	context.move_to(.05, .55)
	context.show_text(gene)
	
	#start legend with title and exon
	context.set_font_size(.06)
	context.move_to(.7, .1)
	context.show_text("Motifs")
	key_pos = .1
	key_pos += .05
	context.set_font_size(.03)
	context.move_to(.7, key_pos)
	context.show_text("EXON")
	
	#add gene position to image
	context.move_to(.05, .6)
	context.show_text(gene_len[gene])
	
	for key in motifs:
		key_count += 1        
		motif_loc_start = .675
		motif_loc_end = .725
		#kmerize gene
		kmerStart = 0
		kmerEnd = kmerStart + len(key)
		context.set_line_width(scale * len(key))
		
		while kmerEnd <= len(gene_info[gene]):  
			#if kmer matches current motif add to image in correct color
			if (re.match(motifs[key], gene_info[gene][kmerStart:kmerEnd])):
				context.move_to((kmerStart * scale) + shift, motif_loc_start)
				context.set_source_rgb(colors[key_count][0], colors[key_count][1], colors[key_count][2])
				context.line_to((kmerStart * scale) + shift, motif_loc_end)
				context.stroke()
				 
			kmerStart += 1
			kmerEnd = kmerStart + len(key)
			
		#input genes into legend
		key_pos += .05
		context.set_source_rgb(colors[key_count][0], colors[key_count][1], colors[key_count][2])
		context.set_font_size(.03)
		context.move_to(.7, key_pos)
		context.show_text(key.upper())
	if args.file_type == ".svg":
		surface.finish()
	elif args.file_type == ".png":
		surface.write_to_png(args.output + gene+ ".png")
		print("png")
	else:
		surface.write_to_png(args.output + gene + args.file_type)