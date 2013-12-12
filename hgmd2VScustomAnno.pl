#!/opt/local/bin/perl -w

use strict;
#use Pod::Usage;  ### doesn't compile with pp on windows.
use Getopt::Long;
use Cwd qw( abs_path );
use File::Copy;

### $VERSION$ ###

=head1 NAME

hgmd2VScustomAnno.pl - convert hgmd data dump into custom annotation

=head1 SYNOPSIS

hgmd2VScustomAnno.pl [--help] [downloaded SNVs and INDELs files]

=head1 DESCRIPTION

Use this script to convert downloaded data from HGMD into a file that
can be loaded into VariantStudio as a custom annotation file. SNVs and
INDELs get downloaded as separate files.

Variants are expressed as you would in a VCF, e.g. REF=GA ALT=G
represents a G deletion.

The script detects the differences in column order and data content
between SNVs and INDELS. Also handles a custom mysql dump

Also compiles a summary list of genes with variants and their
associated diseases.

=head1 AUTHOR

Elliott H. Margulies, Ph.D. -- emargulies@illumina.com
    initiated on 11/01/2013

=head1 OPTIONS

=over 8

=item B<--help>

Display full documentation.

=item B<--vcf>

Generate a VCF file of HGMD variants as well. Requires bgzip.

=item B<--gene-summary>

Generate a gene level summary of HGMD data.

=back

=cut


my $now_time = localtime;
print STDERR "It is now $now_time\n\n";
print STDERR 'This is version $VERSION$' . "\n\n";
print STDERR "Welcome to the converter...\n\n";
print STDERR "This program will convert data downloaded from HGMD into a \n";
print STDERR "custom annotation file that can be imported into " .
    "VariantStudio.\n\n";

my $vcf='';
my $gene_summary='';

GetOptions('vcf' => \$vcf,
	   'gene-summary' => \$gene_summary,);
#	   'help' => sub { pod2usage( verbose => 2) });

my @FILES;
if ($ARGV[0]) {
    @FILES = @ARGV;
} else {
    @FILES = glob("*.txt");
}

if (@FILES) {
    print STDERR "I found the following files to process:\n";
    foreach my $f (@FILES) {
	print STDERR "\t$f\n";
    }
    print STDERR "\n";
} else {
    print STDERR "ERROR: There are no files in this directory for me to " . 
	"process.\n";
    my $path = abs_path();

    print STDERR "ERROR: Please make sure there are files here:\n\t$path\n";
    print STDERR "ERROR: The program is shutting down now without " . 
	"processing anything\n";

    print STDERR "\n  *** Please press the <ENTER> key to exit the program " .
	"and close this window\n";
    <STDIN>;
    print STDERR "bye\n";
    sleep 1;
    exit;
}

### These hashes hold gene summary information, if requested.
my %hgmd_transcript;
my %hgmd_disease;
my %hgmd_omim;
my %hgmd_pmid;

### This hash holds the collected custom annotation data.
my %custom_annotation;


if ($vcf eq '1') {
    open (VCF, "| sort -k1,1 -k2,2n | bgzip -c > hgmd_variants.vcf.gz");
}

foreach my $file (@FILES) {
    print STDERR "working on $file...\n";
    if ($file =~ /\.gz$/) {
	open (IN, "gunzip -c $file |");
    } else {
	open (IN, $file);
    }
 
    while (my $line = <IN>) {
	### skip commented lines
	next if ($line =~ /^\/\// || $line =~ /^\#/ );
	chomp($line);
	my @data = split /\t/, $line;

	### skip header lines that aren't commented
	next if ($data[0] eq 'type' || $data[0] eq 'Type');
	
	### generalize fetching data from multiple sources. Could be
	### either HGMD web downloaded SNVs, Indels, or a custom mysql
	### dump. This is the standard data we need for post
	### processing. get_data() figures out which source we're
	### dealing with and gets the data appropriately.
	my ($class, $ID, $chrom, $start, $stop, $strand, $upstr, 
	    $downstr, $wt, $mut, $gene, $transcript, $HGVSc,
	    $disease, $info, $omim, $pmid) = get_data(@data);

	next if ($chrom eq 'NULL');  ### must have a hg19 coordinate
	next unless ($ID ne '');     ### must have an HGMD ID   

	### INSERTIONS ARE EXPRESSED AS THE BASE BEFORE AND AFTER FOR
	### THE START/STOP RESPECTIVELY
	if ($wt eq 'NULL') {
	    $stop--;
	    $start++;
	}
    
	### express indels as they would in a VCF.
	if (length($wt) > 1 || length($mut) > 1) {
	    make_vcf_allele(\$wt,\$mut,$strand,$upstr,$downstr);
	    $start--;
	}
	
	### clean up output to be all caps
	$wt=uc($wt);
	$mut=uc($mut);

	### convert everything to genomic alleles. coordiantes are
	### already genomic based.
	if ($strand eq '-') {
	    $wt=revcomp($wt);
	    $mut=revcomp($mut);
	}

	$chrom =~ /^chr(\w+)$/;
	my $cnum=$1;

	### prepare VCF output, if requested.
        if ($vcf eq '1') {
	    my $anno = join("|", $ID, $class,$wt,$mut,$gene,$disease,
			    $transcript,$HGVSc,$omim,$pmid,$info);
	    $anno =~ s/NO\ ID//g;
	    $anno =~ s/\ /\_/g;
	    $anno =~ s/\"//g;

	    print VCF join("\t", $cnum, $start, '.', $wt, $mut, '.', '.',
			   "hgmd_variant=$anno") . "\n";
	}
	
	### prepare VariantStudio custom annotation data
	my $cakey = join("\t", "$cnum", $start, $wt, $mut, $ID);
	
	if  ($custom_annotation{$cakey}) {
	    warn "there seems to be a duplicate entry with the same " . 
		"chrom,start,ref,alt,ID\n";
	} else {
	    $custom_annotation{$cakey} = join("\t", $disease, $class,
					      $transcript);
	}
	
	### compile data for gene summary. This is only used if
	### --gene-summary is requested.
	$hgmd_transcript{$gene} .= $transcript . ',';
	$hgmd_disease{$gene}{"$disease:$class"}++;
	if ($omim =~ /^\d+$/) {
	    $hgmd_omim{$gene}{$omim}++;
	}
	if ($pmid =~ /^\d+$/) {
	    $hgmd_pmid{$gene}{$pmid}++;
	}
    }
}

if ($vcf eq '1') {
    close VCF;
}

print_hgmd();

if ($gene_summary eq '1') {
    gene_summary();
}

sub gene_summary {
    open (OUT, ">hgmd_genes.txt");
    foreach my $name (sort {$a cmp $b} keys %hgmd_transcript) {
	my $transcripts = make_unique($hgmd_transcript{$name});
	print OUT $name . "\t" . $transcripts . "\t";
	
	my $line;
	foreach my $disease (keys %{ $hgmd_disease{$name} }) {
	    $line .= $disease . '; ';
	}
	chop($line);
	chop($line);
	print OUT "$line\t";
	
	$line = '';
	foreach my $omim (keys %{ $hgmd_omim{$name} }) {
	    $line .= $omim . '; ';
	}
	chop($line);
	chop($line);
	print OUT "$line\t";
	
	$line = '';
	foreach my $pmid (keys %{ $hgmd_pmid{$name} }) {
	    $line .= $pmid . '; ';
	}
	chop($line);
	chop($line);
	print OUT "$line\n";
    }
    close OUT;
}


sub print_hgmd {
    ### first check that we have something to commit
    unless (%custom_annotation) {
	print STDERR "\nERROR: Your files do not seem to contain any " . 
	    "expected lines.\n";
	print STDERR "ERROR: Please make sure you are providing the " . 
	    "correct files for input\n";
	print STDERR "ERROR: The program is shutting down now without " .
		"making a custom annotation file\n";
	print STDERR "\n *** please press the <ENTER> key to close this " . 
	    "window ***\n";
	<STDIN>;
	sleep 1;
	exit;
    }

    ### print header line
    if ( -e "hgmd_variantstudio_custom_annotation.tsv" ) {
	print STDERR "\nWARN: I found an old file that was processed with " .
	    "this program.\n";

	my $timecode = get_datetime_string();
	
	print STDERR "WARN: Moving old file to:\n";
	print STDERR "\thgmd_variantstudio_custom_annotation_$timecode.tsv\n";
	move("hgmd_variantstudio_custom_annotation.tsv",
	     "hgmd_variantstudio_custom_annotation_$timecode.tsv") 
	    or die "Couldn't move file. Same name? Permissions issue? " . 
	    "Time stamp issue? Please delete the resulting file " .
	    "and try again\n";
    }

    print STDERR "\nPrinting header info to custom annotation file...\n";
    
    open (HGMD, "> hgmd_variantstudio_custom_annotation.tsv");
    print HGMD 'chr';
    foreach my $val (qw(position ref variant annotation annotation2
                        annotation3 annotation4)) {
	print HGMD "\t" . $val;
    }
    print HGMD "\n";
    print STDERR "\nPrinting data to custom annotation file...\n";
    ### print data
    foreach my $cakey (sort {$a cmp $b} keys %custom_annotation) {
	print HGMD $cakey . "\t" . $custom_annotation{$cakey} . "\n";
    }
    
    close HGMD;

    print STDERR "\nFinished printing data\n\nThank you for using this " . 
	"program\n";
    print STDERR "\n *** please press the <ENTER> key to close this " .
	"window ***\n";
    <STDIN>;
    print STDERR "\nbye\n";
    sleep 1;
}



sub get_data {
    ### this subroutine determines the source of the line, so we can
    ### pass it on to another subroutine to correctly parse it.

    my @data = @_;

    ### no null values.
    for my $i (0..23) {
	$data[$i] .= '';
    }
    
    ### initialise variables
    my ($class, $ID, $chrom, $start, $stop, $strand, $upstr, $downstr,
	$wt, $mut, $gene, $transcript, $HGVSc, $disease, $info, $omim, $pmid);

    if ($data[0] =~ /M|S|R|I|D|X|P|N|G|E/) {   ### this is a custome mysql dump
	($class, $ID, $chrom, $start, $stop, $strand, $upstr, $downstr,
	 $wt, $mut, $gene, $transcript, $HGVSc, $disease, $info, $omim, $pmid)
	    = mysql_dump(@data);
    } elsif ($data[0] =~ /4|5|6/) {            ### hgmd indel web dump
	($class, $ID, $chrom, $start, $stop, $strand, $upstr, $downstr,
	 $wt, $mut, $gene, $transcript, $HGVSc, $disease, $info, $omim, $pmid)
	    = web_indel_dump(@data);
    } elsif ($data[0] =~ /1|2|3/) {            ### hgmd snv web dump
        ($class, $ID, $chrom, $start, $stop, $strand, $upstr, $downstr,
         $wt, $mut, $gene, $transcript, $HGVSc, $disease, $info, $omim, $pmid)
	    = web_snv_dump(@data);
    }

    return ($class, $ID, $chrom, $start, $stop, $strand, $upstr, $downstr,
	    $wt, $mut, $gene, $transcript, $HGVSc, $disease, $info, $omim,
	    $pmid);
}

sub web_snv_dump {
    my @data = @_;
    my $class = $data[1];
    my $ID = $data[2];

    my $coord = $data[5];
    $coord =~ /(\w+):(\d+):([+-])/;

    my $chrom = $1;
    my $start = $2;
    my $stop = 'NULL';
    my $strand = $3;

    my $sequence_context = uc($data[11]);

    $sequence_context =~ /^(\w+)\[([ACGT\-]+)\/([ACGT\-]+)\](\w+)$/;

    my $upstr = $1;
    my $downstr = $4;
    my $wt = $2;
    my $mut = $3;

    if ($wt eq '-') {
        $wt = 'NULL';
    }
    if ($mut eq '-') {
        $mut = 'NULL';
    }

    my $gene = $data[8];
    my $HGVSc = $data[6];

    my $transcript = 'NULL';
    if ($HGVSc =~ /^(\w+\.\d+):/) {
	$transcript = $1;
    }

    my $disease = $data[9];
    my $info = $data[12];
    my $omim = 'NULL';
    my $pmid = $data[23];

    return ($class, $ID, $chrom, $start, $stop, $strand, $upstr, $downstr,
            $wt, $mut, $gene, $transcript, $HGVSc, $disease, $info, $omim,
	    $pmid);
}


sub web_indel_dump {
    my @data = @_;

    my $class = $data[1];
    my $ID = $data[19];

    if ($ID eq '') {
	return 'NULL', 'NULL', 'NULL';
    }

    my $coord = $data[6];
    $coord =~ /(\w+):(\d+)/;

    my $chrom = $1;
    my $start = $2;
    my $stop = 'NULL';

    $coord =~ /\(([+-])\)/;

    my $strand = $1;

    my $sequence_context = $data[8];
    
    $sequence_context =~ /^(\w+)\[([ACGT\-]+)\/([ACGT\-]+)\](\w+)$/;

    my $upstr = $1;
    my $downstr = $4;
    my $wt = $2;
    my $mut = $3;

    if ($wt eq '-') {
	$wt = 'NULL';
    }
    if ($mut eq '-') {
        $mut = 'NULL';
    }

    my $gene = $data[2];
    my $HGVSc = $data[4];
    $HGVSc =~ /^(\w+\.\d+):/;

    my $transcript = $1;

    my $disease = $data[3];
    my $info = $data[10];
    my $omim = 'NULL';
    my $pmid = $data[18] . '';

    return ($class, $ID, $chrom, $start, $stop, $strand, $upstr, $downstr,
            $wt, $mut, $gene, $transcript, $HGVSc, $disease, $info, $omim,
	    $pmid);
}


sub mysql_dump {
    my @data = @_;

    my $class = $data[1];
    my $ID = $data[2];
    my $chrom = $data[4];
    my $start = $data[5];
    my $stop = $data[6];
    my $strand = $data[7];
    my $upstr = $data[8];
    my $downstr = $data[9];
    my $wt = $data[10];
    my $mut = $data[11];
    my $gene = $data[12];
    my $transcript = $data[13];
    my $HGVSc = $data[14];
    my $disease = clean_extra_spaces($data[20]);
    my $info = $data[21];
    my $omim = $data[22];
    my $pmid = clean_extra_spaces($data[23]);

    return ($class, $ID, $chrom, $start, $stop, $strand, $upstr, $downstr,
	    $wt, $mut, $gene, $transcript, $HGVSc, $disease, $info, $omim,
	    $pmid);
}


sub clean_extra_spaces {
    my $string = shift;
    chop($string);
    substr $string, 0, 1, "";
    $string =~ s/^\s*//;
    $string =~ s/\s*$//;

    return '"' . $string . '"';
}

sub revcomp {
    my $string = shift;
    $string = reverse($string);
    $string =~ tr/ACGTacgt/TGCAtgca/;
    return($string);
}


sub make_vcf_allele {
    my $ref_wt = shift;
    my $ref_mut = shift;
    my $strand = shift;
    my $upstream = shift;
    my $downstream = shift;

    my $pos = chop($upstream);
    my $neg = substr($downstream,0,1);
    

    if ($$ref_mut eq 'NULL') {
	$$ref_mut = '';
    }
    if ($$ref_wt eq 'NULL') {
        $$ref_wt = '';
    }


    if ($strand eq '+') {
	$$ref_wt = $pos . $$ref_wt;
	$$ref_mut = $pos . $$ref_mut;
    } elsif ($strand eq '-') {
	$$ref_wt = $$ref_wt . $neg;
	$$ref_mut = $$ref_mut . $neg;
    } else {
	die "$strand needs to be + or -. It was neither.\n";
    }
}

sub make_unique {
    my $string = shift;
    unless ($string) {
        return '.';
    }
    my @d = split (",", $string);
    my %u;
    foreach my $val (@d) {
        $u{$val}++;
    }
    my $unique = '';
    foreach my $key (keys %u) {
        $unique .= $key . ',';
    }
    chop($unique);
    return $unique;
}

sub get_datetime_string {
     my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) =
	 localtime(time);
     $sec = sprintf("%02d", $sec);
     $min = sprintf("%02d", $min);
     $hour = sprintf("%02d", $hour);
     my $day = sprintf("%02d", $mday);
     my $month = sprintf("%02d", ($mon+1));
     $year = sprintf("%04d", ($year+1900));
     return join("-", $year, $month, $day) . '_' .
	 join('_', $hour, $min, $sec);
}
