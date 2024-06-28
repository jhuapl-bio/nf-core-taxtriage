#!/usr/bin/env python3

##############################################################################################
# Copyright 2022 The Johns Hopkins University Applied Physics Laboratory LLC
# All rights reserved.
# Permission is hereby granted, free of charge, to any person obtaining a copy of this
# software and associated documentation files (the "Software"), to deal in the Software
# without restriction, including without limitation the rights to use, copy, modify,
# merge, publish, distribute, sublicense, and/or sell copies of the Software, and to
# permit persons to whom the Software is furnished to do so.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
# INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
# PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
# LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
# TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE
# OR OTHER DEALINGS IN THE SOFTWARE.
#

import os
import pandas as pd
import numpy as np
from reportlab.lib import colors
from reportlab.pdfgen import canvas
from reportlab.lib.pagesizes import letter
from distributions import import_distributions, make_vplot, body_site_map
from reportlab.graphics.shapes import Line

from reportlab.platypus import SimpleDocTemplate, Table, TableStyle, Paragraph, Spacer, Image, ListFlowable, ListItem
from reportlab.platypus.flowables import HRFlowable, Flowable

from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import inch

from datetime import datetime
from io import StringIO
from reportlab.lib.colors import Color

import argparse
def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Validate and transform a tabular samplesheet.",
        epilog="Example: python check_samplesheet.py samplesheet.csv samplesheet.valid.csv",
    )
    parser.add_argument(
        "-i",
        "--input",
        metavar="INPUT",
        required=True,
        help="Base pathogen discovery table file, TSV format",
    )
    parser.add_argument(
        "-distributions",
        "--distributions",
        metavar="DISTRIBUTIONS",
        required=False,
        help="TSV file that contains all the distribution information for body sites and organisms",
    )
    parser.add_argument("-a", "--abundance_col", metavar="ABU", required=False, default='% Aligned Reads',
                        help="Name of abundance column, default is abundance")
    parser.add_argument("-x", "--id_col", metavar="IDCOL", required=False, default='Name',
                        help="Name of id column, default is id")
    parser.add_argument("-v", "--version", metavar="VERSION", required=False, default='Local Build',
                        help="What version of TaxTriage is in use")
    parser.add_argument("-s", "--sitecol", metavar="SCOL", required=False, default='Sample Type',
                        help="Name of site column, default is body_site")
    parser.add_argument("-t", "--type", metavar="TYPE", required=False, default='name',
                        help="What type of data is being processed. Options: 'tax_id' or 'name'.",
                        choices=['tax_id', 'name'])
    parser.add_argument(
        "-o",
        "--output",
        metavar="OUTPUT",
        required=True,
        type=str,
        help="Path of output file",
    )

    return parser.parse_args(argv)




# Function to adjust font size based on text length
def adjust_font_size(text, max_length=10, default_font_size=10, min_font_size=6):
    if len(text) > max_length:
        # Calculate new font size (simple linear reduction, could be improved)
        new_size = max(default_font_size - (len(text) - max_length) // 5, min_font_size)
        return f'<font size={new_size}>{text}</font>'
    else:
        return f'<font size={default_font_size}>{text}</font>'

# Ensure cell content that is empty is displayed as a blank space in the PDF
def format_cell_content(cell):
    # Convert NaN or None to an empty string
    if pd.isna(cell):
        return ""
    else:
        # Adjust font size based on content length
        return adjust_font_size(str(cell))



def import_data(inputfile ):
    # Load your TSV data into a DataFrame
    # tsv_data = """
    # Name\tSample\tSample Type\t% Aligned\t% Total Reads\t# Aligned\tIsAnnotated\tPathogenic Sites\tType\tTaxid\tStatus\tGini Coefficient\tMeanBaseQ\tMeanMapQ\tBreadth of Coverage\tDepth of Coverage\tisSpecies\tAnnClass
    # Staphylococcus aureus\tshortreads\tstool\t0.0008\t0.0008\t2\tYes\tstool\tCommensal\t1280\testablished\t0.4\t\t\t\t\tTrue\tNone
    # Klebsiella pneumoniae\tshortreads\tstool\t0.002\t0.002\t20\tYes\t"abscess, stool, skin, urine"\tCommensal\t573\testablished\t0.23\t\t\t\t\tTrue\tDerived
    # Dickeya fangzhongdai\tshortreads\tstool\t0.0002\t0.0002\t2\tYes\t\t\t1778540\tN/A\t0.95\t\t\t\t\tTrue\tDirect
    # Pediococcus acidilactici\tlongreads\tnasal\t0.0005\t0.0005\t5\tNo\t\t\t1254\tN/A\t0.9\t\t\t\t\tTrue\tDerived
    # Neisseria gonorrhoeae\tlongreads\tnasal\t0.0629\t0.0629\t120\tYes\t"blood, oral, stool, urine"\tPathogen\t485\testablished\t0.02\t\t\t\t\tTrue\tDirect
    # Escherichia coli\tshortreads\tstool\t0.01\t0.01\t100\tNo\t\t\t93061\tN/A\t0.48\t\t\t\t\tTrue\tDerived
    # Metabacillus litoralis\tshortreads\tstool\t0.08\t0.08\t800\tNo\t\t\t152268\tN/A\t0.80\t\t\t\t\tTrue\tDirect
    # Fluviibacter phosphoraccumulans\tlongreads\tnasal\t0.0005\t0.0005\t5\tNo\t\t\t1751046\tN/A\t0.96\t\t\t\t\tTrue\tDirect
    # Diaphorobacter ruginosibacter\tlongreads\tnasal\t0.00003\t0.00003\t1\tNo\t\t\t1715720\tN/A\t0.97\t\t\t\t\tTrue\tDerived
    # """.strip()
    # df = pd.read_csv(StringIO(tsv_data), sep='\t')

    # # Simulating additional data
    # np.random.seed(42)
    # # df['Gini Coefficient'] = np.random.uniform(0, 1, df.shape[0])
    # df['MeanBaseQ'] = np.random.uniform(20, 40, df.shape[0])
    # df['MeanMapQ'] = np.random.uniform(30, 60, df.shape[0])
    # df['Breadth of Coverage'] = np.random.uniform(50, 100, df.shape[0])
    # df['Depth of Coverage'] = np.random.uniform(10, 100, df.shape[0])

    df = pd.read_csv(inputfile, sep='\t')

    # sort the dataframe by the Sample THEN the # Reads
    df = df.sort_values(by=[ "Type", "Sample", "# Aligned"], ascending=[False, True, False])
    # trim all of NAme column  of whitespace either side
    df['Name'] = df['Name'].str.strip()
    dictnames = {
        11250: "human respiratory syncytial virus B",
        12814: "human respiratory syncytial virus A",
    }
    # df['Organism'] = df['Name']
    for row_idx, row in df.iterrows():
        if row['Status'] == 'putative':
            #  update index of row to change organism name to bold
            df.at[row_idx, 'Name'] = f'{row.Name}*'
        # change the Name column if the mapnames for taxid is in the dict
        df['Name'] = df[['Name', 'Taxid']].apply(lambda x: dictnames[x['Taxid']] if x['Taxid'] in dictnames else x['Name'], axis=1)
    # replace all NaN with ""
    df = df.fillna("")
    return df

def split_df(df_full):
    # Filter DataFrame for IsAnnotated == 'Yes' and 'No'
    # df_ = df_full[df_full['IsAnnotated'] == 'Yes'].copy()
    df_yes = df_full[~df_full['Type'].isin([ 'Unknown', 'N/A', np.nan, "Commensal", "Opportunistic", "Potential" ] ) ].copy()
    df_opp = df_full[df_full['Type'].isin([ 'Opportunistic', "Potential"])].copy()
    df_comm = df_full[df_full['Type'].isin(['Commensal'])].copy()
    df_unidentified = df_full[(df_full['Type'].isin([ 'Unknown', 'N/A', np.nan, ""] ))].copy()

    df_yes.reset_index(drop=True, inplace=True)
    df_opp.reset_index(drop=True, inplace=True)
    return df_yes, df_opp, df_comm, df_unidentified



# Custom styles for title and subtitle with right alignment
left_align_style = ParagraphStyle(
    name='leftAlign',
    parent=getSampleStyleSheet()['Normal'],
    alignment=0,  # 2 is for right alignment
    # fontSize=12,
    spaceAfter=10,
)
subtitle_style = ParagraphStyle(
    name='leftAlignSubtitle',
    parent=getSampleStyleSheet()['Normal'],
    alignment=0,  # Right align
    # fontSize=12,
    spaceAfter=10,
)

title_style = ParagraphStyle(
    name='leftAlignTitle',
    parent=getSampleStyleSheet()['Title'],
    alignment=0,  # Right align
    # fontSize=14,
    spaceAfter=10,
)
styles = getSampleStyleSheet()
small_font_style = ParagraphStyle(name='SmallFont', parent=styles['Normal'], fontSize=2)
normal_style = styles['Normal']

def prepare_data_with_headers(df, plot_dict, include_headers=True, columns=None):
    data = []
    if not columns:
        columns = df.columns.values[:-1]  # Assuming last column is for plots which should not be included in text headers
    if include_headers:
        headers = [Paragraph('<b>{}</b>'.format(col), styles['Normal']) for col in columns]
        if len(plot_dict.keys()) > 0:
            headers.append(Paragraph('<b>Percentile of Healthy Subject (HHS)</b>', styles['Normal']))  # Plot column header
        data.append(headers)
    for index, row in df.iterrows():
        row_data = [Paragraph(format_cell_content(str(cell)), small_font_style ) for cell in row[columns][:]]  # Exclude plot data
        # row_data = []
        # for col in columns:
        #     cell_content = format_cell_content(str(row[col]))
        #     if col in ["# Aligned to Sample", "Sample"]:
        #         print(col)
        #         # row_data.append(Paragraph(cell_content, small_font_style))
        #         row_data.append(Paragraph(cell_content, small_font_style))
        #     else:
        #         row_data.append(Paragraph(cell_content,small_font_style))
        # Insert the plot image
        if len(plot_dict.keys()) > 0:
            plot_key = (row['Organism'], row['Type'])
            if plot_key in plot_dict:
                plot_image = Image(plot_dict[plot_key])
                plot_image.drawHeight = 0.5 * inch  # Height of the image
                plot_image.drawWidth = 1* inch  # Width of the image, adjusted from your figsize
                row_data.append(plot_image)
        data.append(row_data)

    return data

def return_table_style(df, color_pathogen=False):
    # Start with the basic style
    table_style = TableStyle([
        ('BACKGROUND', (0,0), (-1,0), colors.gray), # header
        ('TEXTCOLOR', (0,0), (-1,0), colors.whitesmoke),
        ('GRID', (0,0), (-1,-1), 1, colors.black),
        ('VALIGN', (0,0), (-1,-1), 'TOP'),

    ])
    if color_pathogen:
        # Placeholder for cells to color (row_index, col_index) format
        cells_to_color = []
        colorindexcol = 2

        sampleindx = df.columns.get_loc('Type')
        # Example post-processing to mark cells
        for row_idx, row in enumerate(df.itertuples(index=False)):
            val = row.Class
            status = row.Status
            sites = row.Locations
            # if nan then set to empty string
            if pd.isna(sites):
                sites = ""
            # Get Sample Type value from row
            sampletype = row[sampleindx]
            if val != "Commensal" and sampletype in sites:
                color = 'lightgreen'
            elif val != "Commensal" and row.AnnClass == 'Derived':
                color = 'lightblue'
            elif val != "Commensal":
                color = "papayawhip"
            elif val == "Commensal" and row.AnnClass == "Derived":
                color = 'lightblue'
            else:
                color = "white"
            # Ensure indices are within the table's dimensions
            style_command = ('BACKGROUND', (colorindexcol, row_idx+1), (colorindexcol, row_idx+1), color)  # Or lightorange based on condition
            table_style.add(*style_command)
    else:
        table_style.add(*('BACKGROUND', (0,1), (-1,-1), colors.white))
    return table_style

def make_table(data, table_style=None):
    # Set table style to the return value of the function

    # Set style
    # Applying the custom style to the title and subtitle
    # Table configuration
    table = Table(data, repeatRows=1)
    # Apply this style to your tables
    table.setStyle(table_style)
    return table
def draw_vertical_line(canvas, doc):
        """
        Draw a vertical line 5% from the left of the page, starting 10% down from the top
        and ending 10% up from the bottom.
        """
        page_width, page_height = letter
        line_x = 0.05 * page_width
        start_y = 0.1 * page_height
        end_y = page_height - (0.1 * page_height)
        canvas.saveState()
        canvas.setStrokeColor(colors.black)
        canvas.setLineWidth(1)
        canvas.line(line_x, start_y, line_x, end_y)
        canvas.restoreState()


def create_report(
    output,
    df_identified,
    df_opportunistic,
    df_unidentified,
    df_commensals,
    plotbuffer,
    version=None
):

    # PDF file setup
    pdf_file = output
    doc = SimpleDocTemplate(pdf_file, pagesize=letter)
    # Placeholder values for version and date
    #### Section to style things up a bit
    # Set the left margin to 10% of the width of a letter size page (8.5 inches)
    left_margin = 0.1 * letter[0]
    right_margin = left_margin / 5
    top_margin = bottom_margin = 0.1 * letter[1]


    # Modify the doc setup to include custom margins
    doc = SimpleDocTemplate(
        pdf_file,
        pagesize=letter,
        leftMargin=left_margin,
        rightMargin=right_margin,
        topMargin=top_margin,
        bottomMargin=bottom_margin
    )
    # version = "1.3.2"  # Example version
    if not version:
        version = "Local Build"
    # get datetime of year-mont-day hour:min
    date = datetime.now().strftime("%Y-%m-%d %H:%M")  # Current date
    # sort df_identified by Alignment Conf
    # filter out so only Class is PAthogen
    df_identified = df_identified.sort_values(by=['Sample', 'Alignment Conf'], ascending=False)
    df_opportunistic = df_opportunistic.sort_values(by=['Sample', 'Alignment Conf'], ascending=False)
    df_identified_paths = df_identified
    df_identified_others = df_commensals
    # df_identified_others = df_identified[df_identified['Class'] != 'Pathogen']
    df_unidentified = df_unidentified.sort_values(by=['Sample', '# Aligned'], ascending=False)
    elements = []
    ##########################################################################################
    ##########################################################################################
    ##### Section to make the Top Table - all annotated commensal or otherwise
    if not df_identified_paths.empty:
        columns_yes = df_identified_paths.columns.values
        # print only rows in df_identified with Gini Coeff above 0.2
        columns_yes = ["Sample (Type)", "Organism", "Class", "# Aligned", "Alignment Conf", "Taxid", "Pathogenic Subsp/Strains"]
        # Now, call prepare_data_with_headers for both tables without manually preparing headers
        data_yes = prepare_data_with_headers(df_identified_paths, plotbuffer, include_headers=True, columns=columns_yes)
        table_style = return_table_style(df_identified_paths, color_pathogen=True)
        table = make_table(
            data_yes,
            table_style=table_style
        )
        # Add the title and subtitle

        title = Paragraph("Organism Discovery Analysis", title_style)
        subtitle = Paragraph(f"This report was generated using TaxTriage <b>{version}</b> on <b>{date}</b> and is derived from an in development spreadsheet of human-host pathogens. It will likely change performance as a result of rapid development practices.", subtitle_style)
        elements.append(title)
        elements.append(subtitle)
        elements.append(Spacer(1, 12))
        elements.append(table)
        elements.append(Spacer(1, 12))  # Space between tables
    # Adding regular text

    styles = getSampleStyleSheet()

    # Adding subtext (you can adjust the style to make it look like subtext)
    subtext_style = styles["BodyText"]
    subtext_style.fontSize = 10  # Smaller font size for subtext
    subtext_style.leading = 12
    subtext_para = Paragraph("Organisms marked with * are putative and have relatively lower references listing their annotations as a pathogen in the given sample types. Classifications of pathogens are described as:", subtext_style)
    elements.append(subtext_para)

    # Create a bullet list
    bullet_list_items = [
        "Primary: Exposure to the agent generally results in a diseased state in both immunocompromised and immunocompetent individuals.",
        "Opportunistic: Exposure to the agent causes a diseased state under certain conditions, including immunocompromised status, wound infections, and nosocomial infections",
        "Commensal: Organisms typically found in the human microbiota.",
        "Potential: Organisms that have been associated with disease states but are not extensively studied",

    ]

    bullet_list = ListFlowable(
        [ListItem(Paragraph(item, subtext_style)) for item in bullet_list_items],
        bulletType='bullet',  # '1' for numbered list
        start='circle'  # 'circle', 'square', etc.
    )

    elements.append(bullet_list)
    # Create an HRFlowable for the horizontal line
    horizontal_line = HRFlowable(width="100%", thickness=1, color=colors.black, spaceBefore=12, spaceAfter=12)

    # Add the horizontal line to the elements
    elements.append(horizontal_line)



    subtext_para = Paragraph("The following information highlights the description for the color combinations for each organism class in the annotated table(s)", subtext_style)
    elements.append(subtext_para)

    bullet_list_items = [
        "Green/White: Direct match for the taxid/organism name with your sample type from the database.",
        "Blue: Derived Pathogenicity from any listed pathogenic strains of a given organism.",
        "Beige: Pathogens annotated in sample type(s) other than your listed one.",
    ]

    # Create a list of bullet items with specified colors
    bullet_colors = [colors.lightgreen, colors.lightblue, colors.papayawhip,  ]
    style = styles['Normal']

    # Create custom ListItems with colored bullets
    custom_list_items = [
        ListItem(Paragraph(item, style), bulletBorder=colors.black,  bulletColor=bullet_colors[idx], )
        for idx, item in enumerate(bullet_list_items)
    ]

    # Create the ListFlowable
    bullet_list = ListFlowable(
        custom_list_items,
        start="square",

        bulletType='bullet'  # '1' for numbered list
    )

    # add horizontal line in reportlab
    elements.append(bullet_list)

    subtext_para = Paragraph("Read amounts are represented as the <b>total number of aligned reads</b> of sufficient mapping quality <b>(% aligned for all reads in sample)</b>", subtext_style)
    elements.append(subtext_para)

    elements.append(horizontal_line)

    ##########################################################################################
    #### Table on opportunistic pathogens
    if not df_opportunistic.empty:
        columns_opp =  ["Sample (Type)", "Organism", "Class", "# Aligned", "Alignment Conf", "Taxid", "Pathogenic Subsp/Strains"]
        data_opp = prepare_data_with_headers(df_opportunistic, plotbuffer, include_headers=True, columns=columns_opp)
        table_style = return_table_style(df_opportunistic, color_pathogen=True)
        table = make_table(
            data_opp,
            table_style=table_style
        )
        # Add the title and subtitle

        Title = Paragraph("Opportunistic Pathogens", title_style)
        elements.append(Title)
        elements.append(Spacer(1, 12))
        elements.append(table)
        elements.append(Spacer(1, 12))  # Space between tables

    ################################################################################################
    ### Table on commensals
    if not df_identified_others.empty:
        columns_yes = df_identified_others.columns.values
        # print only rows in df_identified with Gini Coeff above 0.2
        columns_yes = ["Sample (Type)", "Organism", "Class", "# Aligned", "Alignment Conf", "Taxid"]
        # Now, call prepare_data_with_headers for both tables without manually preparing headers
        data_yes_others = prepare_data_with_headers(df_identified_others, plotbuffer, include_headers=True, columns=columns_yes)
        table_style = return_table_style(df_identified_others, color_pathogen=True)
        table = make_table(
            data_yes_others,
            table_style=table_style
        )
        # Add the title and subtitle
        title = Paragraph("Commensals", title_style)
        subtitle = Paragraph(f"These were identified & were listed as a commensal directly", subtitle_style)
        elements.append(title)
        elements.append(subtitle)
        elements.append(Spacer(1, 12))  # Space between tables

        elements.append(table)
        elements.append(Spacer(1, 12))  # Space between tables
    # # Adding regular text

    elements.append(Spacer(1, 12))
    if not df_unidentified.empty:
        ##########################################################################################
        ### Section to Make the "Unannotated" Table
        second_title = "Unannotated Organisms"
        second_subtitle = "The following table displays the unannotated organisms and their alignment statistics. Be aware that this is the exhaustive list of all organisms (species only) contained within the samples that had atleast one read aligned"
        elements.append(Paragraph(second_title, title_style))
        elements.append(Paragraph(second_subtitle, subtitle_style))

        columns_no = ['Sample', 'Organism','# Aligned', "Alignment Conf" ]
        data_no = prepare_data_with_headers(df_unidentified, plotbuffer, include_headers=True, columns=columns_no)
        table_style = return_table_style(df_unidentified, color_pathogen=False)
        table_no = make_table(
            data_no,
            table_style=table_style
        )
        elements.append(table_no)
    elements.append(Spacer(1, 12))  # Space between tables
    ##########################################################################################
    #####  Build the PDF
    # Adjust the build method to include the draw_vertical_line function
    doc.build(elements, onFirstPage=draw_vertical_line, onLaterPages=draw_vertical_line)

    print(f"PDF generated: {pdf_file}")



def main():
    args = parse_args()
    df_full = import_data(args.input)

    # change column "id" in avbundance_data to "tax_id" if args.type is "name"
    df_full = df_full.rename(columns={args.id_col: args.type})
    df_full = df_full.rename(columns={args.sitecol: 'body_site'})
    df_full = df_full.rename(columns={args.abundance_col: 'abundance'})
    df_full = df_full.dropna(subset=[args.type])
    # df_identified = df_identified[[args.type, 'body_site', 'abundance']]
    # convert all body_site with map
    df_full['body_site'] = df_full['body_site'].map(lambda x: body_site_map(x) )
    # make new column that is # of reads aligned to sample (% reads in sample) string format
    df_full['Quant'] = df_full.apply(lambda x: f"{x['# Aligned']} ({x['abundance']:.2f}%)", axis=1)
    # add body sit to Sample col with ()
    def make_sample(x):
        if not x['body_site']:
            return ""
        else:
            return f"{x['Sample']} ({x['body_site']})"
    df_full['Sample (Type)'] = df_full.apply(lambda x: make_sample(x), axis=1)
    # group on sampletype and get sum of abundance col
    # get the sum of abundance for each sample

    plotbuffer = dict()
    if args.distributions and os.path.exists(args.distributions):
        stats_dict, site_counts = import_distributions(
            args.distributions,
            args.type,
            []
        )

        for index, row in df_full.iterrows():
            # taxid, body_site, stats, args, result_df
            # if taxid and body site not in stats dict then make it empty or 0
            taxidsonly = [key[0] for key in stats_dict.keys()]
            bodysites = [key[1] for key in stats_dict.keys()]
            # if (row['tax_id'], row['body_site']) in dists or row['tax_id'] not in taxidsonly or len(body_sites)== 0:
            if (row[args.type], row['body_site']) not in stats_dict:
                stats = {
                    'mean': 0,
                    'std': 0,
                    "min_abundance":0,
                    "max_abundance": 0,
                    "variance": 0,
                    "body_site": "Unknown",
                    "tax_id": "Unknown",
                    "name": "Unknown",
                    "rank": "Unknown",
                    "abundances": [],
                    "gini_coefficient": 0,
                }
            else:
                stats = stats_dict[(row[args.type], row['body_site'])]
            rank = stats['rank']
            buffer = make_vplot(
                row[args.type],
                stats,
                args.type,
                df_full
            )
            plotbuffer[(row[args.type], row['body_site'])] = buffer
    # convert all locations nan to "Unknown"
    df_full['Pathogenic Sites'] = df_full['Pathogenic Sites'].fillna("Unknown")

    df_full['Alignment Conf'] = df_full['Gini Coefficient'].apply(lambda x: f"{x:.2f}" if not pd.isna(x) else 0)
    print(f"Size of of full list of organisms: {df_full.shape[0]}")
    df_identified, df_opportunistic, df_commensal, df_unidentified= split_df(df_full)
    remap_headers = {
        "Name": "Organism",
        "name": "Organism",
        "# Aligned": "# Reads Aligned to Sample",
        "body_site": "Type",
        "abundance": "% of Aligned",
        "Pathogenic Sites": "Locations",
        "% Reads": "% Reads of Organism",
        "Type": "Class",
        'Quant': "# Aligned",
        "Gini Coefficient": "Gini Coeff",
    }
    df_identified= df_identified.rename(columns=remap_headers)
    df_unidentified= df_unidentified.rename(columns=remap_headers)
    df_commensal = df_commensal.rename(columns=remap_headers)
    df_opportunistic = df_opportunistic.rename(columns=remap_headers)
    version = args.version
    create_report(
        args.output,
        df_identified,
        df_opportunistic,
        df_unidentified,
        df_commensal,
        plotbuffer,
        version,
    )






if __name__ == "__main__":
    main()

