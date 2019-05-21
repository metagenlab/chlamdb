#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-

import os
import shell_command

class Circosviz():


    def __init__(self, fasta_file, bam_file):



        self.working_dir = os.getcwd()
        self.fasta_file = os.path.join(self.working_dir,fasta_file)
        self.bam_file = os.path.join(self.working_dir,bam_file)
        self.prefix = bam_file.split('.')[0]


        self.data_dir = os.path.join(self.working_dir, 'circosviz/data')
        self.etc_dir = os.path.join(self.working_dir, 'circosviz/etc')
        self.sam_file = os.path.join(self.data_dir ,self.prefix+'.sam')

        try:
            os.makedirs(self.data_dir)
        except:
            print '%s already exits' % self.data_dir

        try:
            os.makedirs(self.etc_dir)
        except:
            print '%s already exits' % self.etc_dir

        # convert bam to sam
        print 'converting bam to sam'
        cmd = ' samtools view -h -o %s %s' % (self.sam_file, self.bam_file)
        print cmd
        out, err, code = shell_command.shell_command(cmd)
        if code != 0:
            print out, err

        self.template_config = '''
###################################################################################
#              Basic circos definition file for genome validation                 #
###################################################################################


###################################################################################
############################### Ideogram ##########################################
###################################################################################

karyotype = data/circos.karyotype.txt
chromosomes_units           = 1000
chromosomes_display_default = yes
#chromosomes_order           = 10, 12, 278, 113, 166, 80, 40, 129, 111, 2, 117, 394, 100, 99, 36, 43, 98, 79
#chromosomes_reverse         = 10, 80, 129, 98, 99, 100, 79

<ideogram>

<spacing>
default            = 10u
</spacing>

# thickness and color of ideograms
thickness          = 25p
stroke_thickness   = 1
stroke_color       = black

# the default chromosome color is set here and any value
# defined in the karyotype file overrides it
fill               = yes
fill_color         = black

# fractional radius position of chromosome ideogram within image
radius             = 0.75r
show_label         = no
label_font         = default
label_radius       = dims(ideogram,radius) + 0.175r
label_size         = 30
label_parallel     = no

# show_bands determines whether the outline of cytogenetic bands
# will be seen
show_bands         = yes
band_stroke_thickness = 1

# in order to fill the bands with the color defined in the karyotype
# file you must set fill_bands
fill_bands         = yes
band_transparency  = 1

</ideogram>

###################################################################################
################################# Ticks ###########################################
###################################################################################

show_ticks         = yes
show_tick_labels   = no

<ticks>
tick_label_font    = condensed
radius             = dims(ideogram,radius_outer)
label_offset       = 5p
label_size         = 5p
color              = black
thickness          = 1p

<tick>
spacing           = 100u
size              = 8p
show_label        = no
label_size        = 6p
format            = %d
</tick>
<tick>
spacing           = 10u
size              = 5p
show_label        = no
label_size        = 8p
format            = %d
</tick>

</ticks>

###################################################################################
################################# Plots ###########################################
###################################################################################
<plots>

######################### Coverage plot histogram #################################
<plot>
type               = histogram
thickness          = 1p
color              = black
file               = data/circos.coverage.txt
show		         = yes
r0                 = 1.03r
r1                 = 1.20r
min                = 0
max                = 1000
fill_under         = yes
fill_color         = grey_a5

<rules>
<rule>
condition          = var(value) > 600
color              = red
fill_color         = red
</rule>

<rule>
condition          = var(value) < 20
color              = red
fill_color         = red
</rule>
</rules>

</plot>

######################### GC line plot ############################################
#<plot>
#type		         = line
#show		         = yes
#min               = 30
#max               = 50
#r0                 = 1.43r
#r1                 = 1.47r
#color              = green
#thickness          = 1p
#file               = data/circos.gc.txt
#z                  = 1
# </plot>

######### Coverage plot heatmap - to highlight low coverage areas #################
<plot>
type               = heatmap
file               = data/circos.coverage.txt
show	       	    = yes
r0                 = 1.02r
r1                 = 1.03r
min                = 0
max                = 5000
color              = blues-9-seq
stroke_thickness   = 0
scale_log_base     = 1

<rules>
<rule>
condition          = var(value) > 500
show               = no
</rule>

#<rule>
#condition          = var(value) < 500
#color              = red
#</rule>
</rules>

</plot>

######################## Count of good end connections ############################

<plot>
type               = heatmap
file               = data/circos.count.ends.txt
show	       	    = yes
r0                 = 0.98r
r1                 = 0.99r
#min               = 0
#max               = 2000
color              = green
stroke_thickness   = 0
#scale_log_base    = 1
</plot>

########### Count of connections within same scaffold, but with wrong length ######
<plot>
type               = heatmap
file               = data/circos.count.scontigs.wl.txt
show		         = yes
r0                 = 0.97r
r1                 = 0.98r
#min               = 0
#max               = 5000
color              = blues-9-seq
stroke_thickness   = 0
#scale_log_base    = 1

<rules>
<rule>
condition          = var(value) > 50
color              = red
</rule>

<rule>
condition          = var(value) < 10
show               = no
</rule>
</rules>

</plot>

################## Count of non-end connections between scaffolds #################
<plot>
type               = heatmap
file               = data/circos.count.dcontigs.txt
show		         = yes
r0                 = 0.96r
r1                 = 0.97r
#min               = 0
#max               = 5000
color              = blues-9-seq
stroke_thickness   = 0
#scale_log_base    = 1

<rules>
<rule>
condition          = var(value) > 10
color              = blue
</rule>

<rule>
condition          = var(value) < 50
show               = no
</rule>
</rules>

</plot>

###################################################################################
############################## FRCbam stats #######################################
###################################################################################

#Each parameter from FRCbam is assigned a number from 1 to 14.
#Use these numbers to select which FRCbam parameters to show (use rules).

#STRECH_PE          = 1
#COMPR_PE           = 2
#LOW_NORM_COV_PE    = 3
#LOW_COV_PE         = 4
#HIGH_COV_PE        = 5
#HIGH_NORM_COV_PE   = 6
#HIGH_OUTIE_PE      = 7
#HIGH_SINGLE_PE     = 8
#HIGH_SPAN_PE       = 9
#COMPR_MP           = 10
#HIGH_OUTIE_MP      = 11
#HIGH_SINGLE_MP     = 12
#HIGH_SPAN_MP       = 13
#STRECH_MP          = 14

########################### Show all FRCbam features###############################
#<plot>
#type               = heatmap
#file               = data/frc.tracks.txt
#show		         = yes
#r0                 = 1.01r
#r1                 = 1.015r
#color              = red
#stroke_thickness   = 0

#</plot>

############################### COMPR_PE  ########################################
#<plot>
#type               = heatmap
#file               = data/frc.tracks.txt
#show		         = yes
#r0                 = 1.015r
#r1                 = 1.02r
#color              = blue
#stroke_thickness   = 0

#<rules>
#<rule>
#condition          = var(value) != 2   #Change which FRCbam parameter is shown here
#show               = no
#</rule>
#</rules>

#</plot>

</plots>


###################################################################################
################################# Links ###########################################
###################################################################################
<links>
bezier_radius      = 0.1r
radius             = 0.95r

####################### Links of good end connections #############################
<link chain1>
file               = data/circos.ends.txt
thickness          = 1
color              = black_a5
show               = yes

</link chain1>

### Links of connections within same scaffold, but with wrong length ##############
<link chain2>
file               = data/circos.scontigs.wl.txt
thickness          = 1
color              = red_a5
show               = no

</link chain2>

############### Links of non-end connections between scaffolds ####################
<link chain3>
file               = data/circos.dcontigs.txt
thickness          = 1
color              = blue_a5
show               = yes
</link chain3>

</links>

###################################################################################
#################### Additional system settings ###################################
###################################################################################

<image>
<<include etc/image.conf>>
</image>

# includes etc/colors.conf
#          etc/fonts.conf
#          etc/patterns.conf
<<include etc/colors_fonts_patterns.conf>>

# system and debug settings
<<include etc/housekeeping.conf>>

anti_aliasing*     = no

        '''


    def generate_circos_data_files(self):
        os.chdir(self.data_dir)
        cmd = 'circosviz.pl -i %s -f %s -e 500 -m 3000 -a 125 -b 1000 -p 1000' % (self.sam_file, self.fasta_file)
        out, err, code = shell_command.shell_command(cmd)
        if code != 0:
            print out, err
        os.chdir(self.working_dir)

    def execute_circos(self):

        os.chdir(self.etc_dir)

        self.config_path = 'circosviz.config'
        self.plot_path = 'circosviz_plot.svg'

        with open(self.config_path, 'w') as f:
            f.write(self.template_config)

        cmd = 'circos -noparanoid -conf %s -outputfile %s' % (self.config_path, self.plot_path)
        out, err, code = shell_command.shell_command(cmd)
        if code != 0:
            print out, err

        os.chdir(self.working_dir)

if __name__ == '__main__':
    import argparse
    from Bio import SeqIO
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", '--input_fasta', type=str, help="input fasta file")
    parser.add_argument("-b", '--input_bam', type=str, help="input bam file")

    args = parser.parse_args()
    viz = Circosviz(args.input_fasta, args.input_bam)
    print 'running circosviz'
    viz.generate_circos_data_files()
    print 'running circos'
    viz.execute_circos()

