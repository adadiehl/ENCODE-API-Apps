#!/usr/bin/perl

# search_human-mouse.pl: An extension of search_encode.pl specialized to locate
# and download files for which data is available in human and mouse (or two
# arbitrary other species -- see options). Tabular metadata is produced by
# default and json for each file is stored optionally.
#
# Based on:
#
# search_encode.pl: A general framework to query and retrieve data from the
# ENCODE repository using the REST API. Copyright (C) 2015, Adam Diehl,
# adadiehl@umich.edu
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

use strict;
use warnings;
use JSON;
use WWW::Mechanize;
use Getopt::Long;
use Encode qw(encode_utf8 decode_utf8);
use Digest::MD5 qw(md5_hex);

my $usage =
"\nsearch_human-mouse.pl

A framework to query and retrieve data from the ENCODE repository for datasets
shared between human and mouse. Based on search_encode.pl. Tabular metadata
for search results is always returned (one line for each file found). Direct 
data download and saving of json data are available as options. Currently only
allows comparison by target of assay.

USAGE:

search_human-mouse.pl <search term 1,...,search term N> [OPTIONS]

<search term 1,...,search term N>
    Comma-delimited string of search terms. Can be cell line, tissue type, etc.
    Try to use as precise a term as possible for best results!

OPTIONS:

--species <scientific name 1,...,scientific name N>
    Comma-delimited list of scientific names for all query terms, in
    the same order. Each must be valid scientific name of a species
    present in the ENCODE repository. See ___ for a list. Default:
    \"Homo sapiens,Mus musculus\".

--params
    Used to supply additional search parameters to narrow the results.
    Supplied as key=value pairs corresponding to columns in the
    given database table. See the schema at the address above. Some columns
    have sub-keys that can be searched. See the supplemental table SX at
    http://address_for_supplement.com/supplement.xls for a listing.

    Possible keys of interest for experiments:

        biosample_term_name: tissue/cell used in assay
        assay_term_name: type of assay (e.g., ChIP-seq)
        target.investigated_as: type of target (e.g, transcription factor)
        target.label: name of factor probed (e.g., CTCF)
        status: status of experiment (released/revoked)

--download
    Download files associated with the results. Type(s) of files may be
    specified with the --file-type flag. Currently only tested for results
    of type \"experiment\". May behave unpredictably for other types!

--output-type <type_1,...,type_n>
    Used with --download, limits results to files of given type(s).
    Multiple types may be supplied as a comma-separated list. File types are
    matched against the \"output_type\" column of the file records, and
    available values will vary depending on the experiment. Common values are:

    * reads: raw sequnence reads (usually fastq format)
    * alignments: sequence alignment (usually BAM format)
    * peaks: processed peaks from some peak caller (often bed or bigBed)
    * signal: raw signal (often bigWig format)

--file-format <format_1,...,format_n>
    Used with --output-type, restrict downloads to a given file format, e.g.,
    bam, bigBed. Available values vary depending on the specific experiment.
    This MUST contain a value for every file type supplied with --output-type
    or it will throw an error. Specify \"0\" for all types not requiring a
    file type.

--file-format-type <type_1,...,type_n>
    Used with --output-type, restrict downloads to a given file format type,
    e.g., \"broadPeak\" or \"narrowPeak\". See notes for --file-format.

--download-path <path>
    Location (absolute or relative to workding directory) to write downloaded
    files.

--out-root <root>
    String to prepend to output files.

--save-json
    Save the JSON metadata for every file downloaded. This can be useful for
    identifying individual files and properties later on when retrieving
    multiple files, especially if further programmatic processing is
    anticipated.

--help
    Display this message


NOTES:

* A search string of \"*\" will return all results matching the rest of the
parameters. This can be useful if you want to retrieve the same type of data
for all available cell types, for instance.

* Most table columns in the ENCODE database tables are visible to the search
algorithm and you may chain as many as you would like onto your query.
Please consult the schema, at https://www.encodeproject.org/profiles/graph.svg
to see the columns for each data table.

* Columns in tables referenced by the primary table may also be searched. (All
links are described in the schema) For instance, to search for experiments for
human only based on scientific name, you need to use the \"replicates\"
coumn, which (indirectly) references the \"organism\" table:
&replicates.library.biosample.donor.organism.scientific_name=Homo sapiens

* All column names and values are CASE-SENSITIVE!

* Database values may include white-space. To include a term with white space
in your query, put it in quotes when supplying it at the command line.

* The same key may be used multiple times with different values and results
matching any of the values will be returned. This can be useful to retrieve
data for multiple cell types or transcription factors, for example.

* An excellent tutorial on using the API and navigating the JSON data is
available at: https://www.encodeproject.org/help/rest-api/


EXAMPLE:

Find all transcription-factor ChIP-seq datasets for factors with data in both
human K562 and mouse MEL cells, and download the bigBed peak annotation files
for all of them.

./search_human-mouse.pl K562,MEL --params \"&assay_term_name=ChIP-seq&target.investigated_as=transcription factor\" --out-root chipseq --download --output-type peaks --file-format bigBed


CREDITS AND LICENSE:

Copyright (C) 2017, Adam Diehl

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

Any questions/comments regarding this program may be directed to Adam Diehl:
adadiehl\@umich.edu

Originally published in \"Deciphering ENCODE\". If you use this tool in your 
research, please cite:

Diehl AD and Boyle AP (2016) Deciphering ENCODE.
Trends in Genetics. Volume 32, Issue 4, 238-249.

\n";

if ($#ARGV + 1 < 2) {
    die "$usage\n";
}

my @search_terms = split /,/, $ARGV[0];

my $help = 0;
my $download = 0;
my $output_type_str;
my $file_format_str;
my $file_format_type_str;
my $download_path = "";
my $out_root = "";
my $save_json = 0;
my $species_str = "Homo sapiens,Mus musculus";
my $params = "";

GetOptions (
    "help" => \$help,
    "download" => \$download,
    "output-type=s" => \$output_type_str,
    "file-format=s" => \$file_format_str,
    "file-format-type=s" => \$file_format_type_str,
    "download-path=s" => \$download_path,
    "out-root=s" => \$out_root,
    "save-json" => \$save_json,
    "species=s" => \$species_str,
    "params=s" => \$params
    );


# Check for help option and exit with usage message if found.
if ($help) {
    die "$usage\n";
}

# Check for other options and do some sanity checks...
if (defined($output_type_str) && !$download) {
    die "\n--output-type only works with --download. Try --help.\n";
}
if (defined($file_format_type_str) && !defined($file_format_str)) {
    die "\n--file-format-type requires --file-format. Try --help.\n";
}
if ($download_path ne "" && !$download) {
    die "\n--download-path requires --download. Try --help\n";
}
if ($download && !defined($output_type_str)) {
    print STDERR "\nWarning: --download will retrieve all accepted files associated with search results. If this is not what you want, use --file-type to specify which type(s) you would like to download. See --help.\n\n";
}

# Check the download path for a trailing slash and add one if needed
if ($download_path ne "") {
    if ($download_path !~ m/\/$/) {
	$download_path .= '/';
    }
}

# If file types have been specified, process the string(s).
my @output_types;
if (defined($output_type_str)) {
    @output_types = split /,/, $output_type_str;
}

my @file_formats;
if (defined($file_format_str)) {
    @file_formats = split /,/, $file_format_str;
    if ($#file_formats != $#output_types) {
        die "\n--file-format and --output-type must have the same number of fields!\n";
    }
}

my @file_format_types;
if (defined($file_format_type_str)) {
    @file_format_types = split /,/, $file_format_type_str;
    if ($#file_format_types != $#output_types) {
        die "\n--file-format-type and --output-type must have the same number of fields!\n";
    }    
}

# Set up the virtual browser
my $mech = WWW::Mechanize->new(
    ssl_opts => {
        verify_hostname => 0,
        # Quick and dirty hack to avoid errors about missing SSL certificates
    },
    );

# Prepare the species array
my @species = split /,/, $species_str;

#################################
# Step 1: Run the primary queries for each species against the ENCODE Portal

# Get search results for each species/term combination and store in an array
# of hashes.
my @data;
for (my $i = 0; $i <= $#search_terms; $i++) {
    push @data, &get_search_data($mech, $species[$i], $search_terms[$i],
				 $params);
}


#################################
# Step 2: Find the intersection between the sets of experiments

# Cycle through the first two experiments and store all data for
# which there is a matching key in both. Use the stored intersection
# as the basis for the next comparison, and so on...
print STDERR "Finding intersecting terms...\n";

my %int = compare_res_hashes($data[0], $data[1]);
for (my $i = 2; $i <= $#data; $i++) {
    %int = compare_res_hashes(\%int, $data[$i]);
}

my $n_factors = keys(%int);
print STDERR "Found $n_factors terms with data in both species.\n";
my @intersection = stuff_array(\%int);

#################################
# Step 3: Find the files we want within experiments in the intersection

# Cycle through @intersection, build the metadata table and, if called for,
# create the list to download data files and/or json data.
printf STDERR "Retrieving data for %d experiments...\n", ($#intersection + 1);
my @downloads;
my @metadata;
my $n_files = 0;
my $n_exp = $#intersection + 1;
foreach my $row (@intersection) {

    my %result = %{$row};

    if ($download) {
	
	my @files = @{$result{files}};
	foreach my $file (@files) {
	    my $url = 'http://www.encodeproject.org' . $file . '?format=json';
	    
	    # Get the JSON data
	    my $file_json = &get_json($mech, $url);
	    
	    # Examine the JSON to see if the file matches our criteria
	    my $use_rec = 1; # Download all files by default                    
	    for (my $i = 0; $i <= $#output_types ; $i++) {
		
		# If output_type does not match, reject the file
		if (@output_types) {
		    if (${$file_json}{output_type} ne $output_types[$i]) {
			$use_rec = 0;
		    }
		}
		
		# If file_format does not match, reject the file
		if (@file_formats) {
		    if (${$file_json}{file_format} ne $file_formats[$i]) {
			$use_rec = 0;
		    }
		}
		
		# If file_format_type does not match, reject the file
		if (@file_format_types) {
		    if (${$file_json}{file_format_type}
			ne $file_format_types[$i]) {
			$use_rec = 0;
		    }
		}
	    }
	    
	    if ($use_rec) {
		push @downloads, $file_json;
		$n_files++;

		# Produce a file-centric metadata row
		my $controls_str = join ",", &nopaths($result{possible_controls});
		my $documents_str = join ",", &nopaths($result{documents});
		
		my @row = (&nopath(${$file_json}{href}),
			   ${$file_json}{accession}, ${$file_json}{output_type},
			   ${${$file_json}{replicate}}{biological_replicate_number},
			   ${${$file_json}{replicate}}{technical_replicate_number},
			   $result{dataset_type}, $result{biosample_term_name},
			   $result{biosample_type}, $result{assay_term_name},
			   &nopath($result{target}), $result{status},
			   $result{date_released}, &nopath($result{lab}),
			   $result{accession}, $controls_str, $documents_str);
		push @metadata, \@row;
	    }
	}
    } else {
	# Produce a experiment-centric metadata row
	my $files_str = join ",", &nopaths($result{files});
	my $controls_str = join ",", &nopaths($result{possible_controls});
	my $documents_str = join ",", &nopaths($result{documents});
	
	my @row = ($result{accession}, $result{dataset_type},
		   $result{biosample_term_name}, $result{biosample_type},
		   $result{assay_term_name}, &nopath($result{target}),
		   $result{status}, $result{date_released},
		   &nopath($result{lab}), $files_str, $controls_str,
		   $documents_str);
	push @metadata, \@row;
    }	    
}


######################################
# Step 4: Download the files and metadata

# Set up the metadata file handle and header
my $DATA;
my $outfile = $download_path . $out_root . '.metadata';
# Make sure file name does not start with a . (i.e., no path or root
# were supplied) 
$outfile =~ s/^\.//;
open $DATA, '>', $outfile;

# Build the header
my @header;
if ($download) {
    # File-centric header
    @header = ("file", "accession", "output_type", "biological_replicate",
               "technical_replicate", "dataset_type", "biosample_term_name",
               "biosample_type", "assay_term_name", "target", "status",
               "date_released", "lab", "parent_accession", "possible_controls",
               "documents");
} else {
    # Experiment-centric header
    @header = ("accession", "dataset_type", "biosample_term_name",
               "biosample_type", "assay_term_name", "target", "status",
               "date_released", "lab", "files", "possible_controls",
               "documents");
}
&print_array(\@header, "\t", $DATA);

# Download the files
foreach my $file_json (@downloads) {
    
    &download_file($mech, $file_json, $download_path, $out_root);
    if ($save_json) {
        &save_json($file_json, $download_path, $out_root);
    }
}

# Print the metadata
foreach my $row (@metadata) {
    # Check for any undefined values
    for (my $i = 0; $i <= $#{$row}; $i++) {
        if (!defined(${$row}[$i])) {
            ${$row}[$i] = '.';
        }
    }
    &print_array($row, "\t", $DATA);
}

close $DATA;

print STDERR "Retrieved $n_files files for $n_exp experiments.\n";
print STDERR "Metadata written to $outfile.\nDone\n";


##########################################
# Subroutines

sub get_search_data {
    my $mech = $_[0];
    my $species = $_[1];
    my $search_str = $_[2];
    my $params = $_[3];

    # Build the query URL from the command line args and parameters
    my $URL = 'http://www.encodeproject.org/search/?searchTerm=';
    $URL .= $search_str;
    $URL .= '&replicates.library.biosample.donor.organism.scientific_name=';
    $URL .= $species;
    $URL .= "&type=experiment";
    $URL .= $params;

    # &limit=all: Return all results, not just the first 25
    # &frame=object: Return all non-null fields within the JSON data
    # &format=json: Return search data in JSON format
    $URL .= '&limit=all&frame=object&format=json';

    print STDERR "Query URL: $URL\n";

    # Get the JSON from the URL
    my $json = &get_json($mech, $URL);

    # Print status and result count to the terminal
    print STDERR "${$json}{notification}: ${$json}{total} results found for $species.\n\n";
    
    if (${$json}{total} == 0) {
        # No results found. Return -1
	return -1;
    }

    # Hash to contain results.
    my %results;

    # Individual search results are within the embedded array @graph.
    # Cycle through the results and store in the results hash.
    foreach my $row (@{${$json}{'@graph'}}) {
	# Each row of the array contains a hash reference to a search result.
	# Make a copy of the hash for convenience.
	my %result = %{$row};
        # Get the target string, without the path.
	my $target = &nopath($result{target});
        # Strip the species suffix from the target. The goal here is to make
        # target names comparable between species.
	$target =~ s/-\S+//;
#	print STDERR "$target\n";
	
        # Store a reference to the JSON as an entry to the array wihtin the
	# results hash.
	push @{$results{$target}}, $row;
    }
    
    return \%results;
}

sub get_json {
    # Retrieve JSON data from a URL and return it as a data structure.
    my ($mech, $url) = @_;
    print STDERR "$url\n";
    $mech->get($url);    
    my $json = decode_json($mech->content);
    return $json;
}

sub save_json {
    # Save the JSON data for a file we have downloaded.
    my $json = $_[0];
    my $download_path = $_[1];
    my $out_root = $_[2];

    my $outfile = $download_path . $out_root
	. '.' . ${$json}{accession} . '.json';
    open my $OUT, '>', $outfile;
    print $OUT to_json($json, {pretty=>1});
    close $OUT;
}

sub download_file {
    # Given a file record, download the file from the ENCODE repository and
    # verify the MD5 checksum.
    my $mech = $_[0];
    my $json = $_[1];
    my $download_path = $_[2];
    my $out_root = $_[3];

    my $url = "https://www.encodeproject.org" . ${$json}{href};
    $mech->get( $url );

    print STDERR "Found a matching record at $url. Retrieving data...\n";
    
    # Make sure we got the file we were expecting by verifying the md5 checksum
    print STDERR "\tVerifying Checksum...\n";
    my $checksum = md5_hex($mech->content);
    if ($checksum ne ${$json}{md5sum}) {
	print STDERR "\tERROR: Checksums do not match. Skipping file!\n";
    } else {
	# Checksums match. Save file to the specified
	# path and file name.
	my $outfile = $download_path . $out_root
	    . '.' . ${$json}{accession} . '.' .
	    ${$json}{file_format};
	print STDERR "\tSaving file to $outfile...\n";
	$mech->save_content($outfile);
    }    
}

sub print_array {
    # A generic function to print an array to a file handle as a string with
    # the supplied delimiting.
    my $array = $_[0];
    my $delim = $_[1];
    my $fh = $_[2];

    my $out_str = join $delim, @{$array};
    print $fh "$out_str\n";
}

sub nopaths {
    # A generic function to strip the path from file or directory names stored
    # in an array.
    my @array = @{$_[0]};
    
    my @out;
    if ($#array < 0) {
	push @out, ".";
    } else {
	foreach my $str (@array) {
	    my $s = &nopath($str);
	    push @out, $s;
	}
    }
    return @out;
}

sub nopath {
    # Strip the path from a file or directory name.
    my $str = $_[0];
    if (!defined($str)) {
        return ".";
    }
    $str =~ s/\/$//;
    $str =~ s/\/\S+\///;
    return $str;
}

sub compare_res_hashes {
    my ($res_1, $res_2) = @_;
    my %intersection;
    foreach my $key (sort(keys(%{$res_1}))) {
	#    print STDERR "$key\n";
	if (exists($res_2->{$key})) {
	    $intersection{$key} = [@{$res_1->{$key}}, @{$res_2->{$key}}]
	}
    }
    return %intersection;
}

sub stuff_array {
    my ($int) = @_;

    my @res;
    foreach my $key (sort(keys(%{$int}))) {
	push @res, @{$int->{$key}};
    }
    return @res;
}
