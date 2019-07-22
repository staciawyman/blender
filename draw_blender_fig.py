import svgwrite
import sys
import os
from operator import itemgetter

# Stacia Wyman 21 July 2019
# Python script to draw alignments of BLENDER hit to the target sequence 
# Adapted from visualize.py from https://github.com/aryeelab/guideseq


boxWidth = 10
box_size = 16

colors = {'G': '#F5F500', 'A': '#FF5454', 'T': '#00D118', 'C': '#26A8FF', 'N': '#B3B3B3'}


# Read in filtered off-target list, grab relevant columns, sort by disco score, then mismatches, 
# then genomic coordinate
def parseSitesFile(infile):

    offtargets = []
    with open(infile, 'r') as f:
        f.readline()
        for line in f:
            line_items = line.split('\t')
            # Skip header line
            if "Discoscore" in line:
                continue
            # concatenate PAM to putative guide sequence
            offtarget_sequence = line_items[6]+line_items[5]
            disco_score = line_items[2] 
            coords = line_items[0]
            mismatches = line_items[7]
            if offtarget_sequence != '':
                offtargets.append({'seq': offtarget_sequence.strip(),
                                   'disco_score': int(disco_score.strip()),
                                   'coords': coords.strip(),
                                   'mismatches': int(mismatches.strip())} )
    # Sort by disco score, then by mismatches
    offtargets = sorted(offtargets, key=lambda x: (-x['disco_score'], x['mismatches'], x['coords']))
    return offtargets 


def visualizeOfftargets(infile, outfile, target_guide, title=None):

    output_folder = os.path.dirname(outfile)
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # Get offtargets array from file
    offtargets  = parseSitesFile(infile)

    # Initiate canvas
    dwg = svgwrite.Drawing(outfile + '.svg', size=('8.5in','11in'))
    #dwg = svgwrite.Drawing(outfile + '.svg', profile='full')

    # Define top and left margins
    x_offset = 20
    y_offset = 50

	# Header line for figure, use title if there is one, otherwise...bam file name?
    if title is None:
        dwg.add(dwg.text("DISCOVER-Seq Off Targets", insert=(x_offset, 30), style="font-size:20px; font-family:Arial"))
    else:
        # Define top and left margins
        dwg.add(dwg.text("DISCOVER-Seq Off Targets: "+title, insert=(x_offset, 30), style="font-size:20px; font-family:Arial"))

    x = x_offset+box_size*(len(target_guide)-3)+7
    y = y_offset - 2
    dwg.add(dwg.text("PAM", insert=(x,y), style="font-size:14px; font-family:Courier"))
    # Draw reference sequence row
    for i, c in enumerate(target_guide):
        y = y_offset
        x = x_offset + i * box_size
        dwg.add(dwg.rect((x, y), (box_size, box_size), fill=colors[c]))
        dwg.add(dwg.text(c, insert=(x + 3, y + box_size - 3), fill='black', style="font-size:15px; font-family:Courier"))

    y = y_offset + box_size + 5
    dwg.add(dwg.text("DISCO",insert=(x_offset+box_size*(len(target_guide)+.5),y - 15),style="font-size:15px; font-family:Arial"))
    dwg.add(dwg.text("Score", insert=(x_offset+box_size * (len(target_guide) + .5), y), style="font-size:15px; font-family:Arial"))
    dwg.add(dwg.text("Mismatches", insert=(x_offset+box_size * (len(target_guide) + 4.25), y), style="font-size:15px; font-family:Arial"))
    dwg.add(dwg.text("Coordinates", insert=(x_offset+box_size * (len(target_guide) + 10.7), y), style="font-size:15px; font-family:Arial"))

    # Draw aligned sequence rows
    y_offset += 10  # leave gap after target sequence
    for j, seq in enumerate(offtargets):
        y = y_offset + j * box_size
        for i, c in enumerate(seq['seq']):
            x = x_offset + i * box_size
            if c == target_guide[i] or target_guide[i] == 'N':
                dwg.add(dwg.text(u"\u2022", insert=(x + 4.5, 2 * box_size + y - 4), fill='black', style="font-size:10px; font-family:Courier"))
            else:
                dwg.add(dwg.rect((x, box_size + y), (box_size, box_size), fill=colors[c]))
                dwg.add(dwg.text(c, insert=(x + 3, 2 * box_size + y - 3), fill='black', style="font-size:15px; font-family:Courier"))

        y = y_offset + box_size * (j + 2) - 2
        sa = "font-size:14px; font-family:Arial"
        sc = "font-size:14px; font-family:courier"
        # If this is the target sequence, output it in red
        if int(seq['mismatches']) == 0:
            mismatch_text = dwg.text("Target", insert=(box_size * (len(target_guide) + 6.25) , y), fill='red', style=sa)
            disco_score_text = dwg.text(str(seq['disco_score']), insert=(box_size * (len(target_guide) + 2) , y), fill='red', style=sa)
            coords_text = dwg.text(str(seq['coords']), insert=(box_size * (len(target_guide) + 12) , y), fill='red', style=sc)
        else:
            mismatch_text = dwg.text(seq['mismatches'], insert=(box_size * (len(target_guide) + 7) , y), fill='black', style=sa)
            disco_score_text = dwg.text(str(seq['disco_score']), insert=(box_size * (len(target_guide) + 2) , y), fill='black', style=sa)
            coords_text = dwg.text(str(seq['coords']), insert=(box_size * (len(target_guide) + 12) , y), fill='black', style=sc)
        dwg.add(disco_score_text)
        dwg.add(mismatch_text)
        dwg.add(coords_text)

    dwg.save()


def main():
    if len(sys.argv) >= 4:
        if len(sys.argv) == 5:
            title = sys.argv[4]
        else:
            title = None
        visualizeOfftargets(sys.argv[1], sys.argv[2], sys.argv[3], title=title)
    else:
        print ('Usage: python draw_blender_fig.py INFILE OUTFILE GUIDE [TITLE]')


if __name__ == '__main__':
    main()
