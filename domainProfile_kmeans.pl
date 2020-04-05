#!/usr/bin/env perl 
#===============================================================================
#
#         FILE:  domainProfile_kmeans.pl
#
#        USAGE:  ./domainProfile_kmeans.pl  
#
#       AUTHOR:  Daniel Lang (dl), 
#      VERSION:  1.0
#===============================================================================
use strict;
use warnings;
# AUTHOR:   lang
# CREATED:  23:22:00 16/02/2012
# MODIFIED: 23:22:00 16/02/2012

use FindBin;
use lib $FindBin::Bin ."/lib";
use Getopt::Long;
use Bio::Graphics;
use Bio::SeqIO;
use Bio::Graphics::Feature;
use GD::SVG;
use Algorithm::NeedlemanWunsch;
use Statistics::R;
use Bio::Matrix::Generic;
use File::Temp qw/tempfile/;

my $USAGE   = "$0 -h | -o outpath -n profile_name [-l] [-p] [-k kmeans_cluster_with_max_k -c criterion_for_k_selection[ch|asw]] [-r ratio_considered_as_representative -q qvalue_cutoff_for_outliers -a qvalue_cutoff_for_absent]  -s seq_file [-f format(fasta)] [-i ignore_motif_regex] JALVIEW_ANNOT_FILE\n"; 

my ($opts_h, $opts_l, $opts_o, $opts_n, $opts_p, $opts_f, $opts_s, $opts_k, $opts_r, $opts_q, $opts_a, @opts_i,$criterion);
GetOptions(	'help'			=> \$opts_h,
			'label'			=> \$opts_l,
			'outpath=s'		=> \$opts_o,
			'name=s'		=> \$opts_n,
			'pairs'			=> \$opts_p,
			'format=s'		=> \$opts_f,
			'seqfile=s'		=> \$opts_s,
			'kmeans=i'		=> \$opts_k,
			'ratio=f'		=> \$opts_r,
			'qvalue=f'		=> \$opts_q,
			'absent=f'		=> \$opts_a,
			'ignore=s'		=> \@opts_i,
			'criterion=s'	=> \$criterion,
);

if ($opts_o && ! -d $opts_o) {
    mkdir($opts_o) or die $!;
}

$opts_o="." unless $opts_o && -d $opts_o;

my $r;
my ($P,$PN);
if ($opts_p) {
    open($P,">$opts_o/$opts_n.pair_scores.txt") or die "Cannot create $opts_o/$opts_n.pair_scores.txt! $!\n";
    open(OUT,">$opts_o/$opts_n.clusters.pair_scores.txt") or die "Cannot create $opts_o/$opts_n.clusters.pair_scores.txt! $!\n";
	$PN			= "$opts_o/$opts_n.pair_scores.txt";
}
else {
	($P,$PN) 	= tempfile();
}
if ($opts_r) {
    open(REPORT,">$opts_o/$opts_n.cluster_outlier_domains.txt") or die "Cannot create $opts_o/$opts_n.cluster_outlier_domains.txt! $!\n";
    $r  = Statistics::R->new();
    $opts_r = $opts_r >= 0 && $opts_r<= 1 ? $opts_r : 0.9;
    $opts_q = $opts_q >= 0 && $opts_q<= 1 ? $opts_q : 0.05;
    $opts_a = $opts_a >= 0 && $opts_a<= 1 ? $opts_a : 0.9;
}
if ($criterion) {
	die $USAGE unless $criterion =~ /^asw|ch$/;
}
else {
	$criterion="ch";
}

die $USAGE if $opts_h or !(@ARGV && $opts_s && -f $opts_s && $opts_n);

my (%motifs,%k,%L,@ord,%dom_freq, %ML);

my $gap_open 		= -2;
my $gap_extension 	= -1 ;
my $comma 			= 3;

#############READ SEQFILE
my $in = Bio::SeqIO->new(-file=>$opts_s,-format=>$opts_f ? $opts_f : "fasta");

while (my $seq=$in->next_seq) {
    my $id  = $seq->display_id;
    $id     =~ s/\/\d+-\d+$//; # remove trailing alignment coords
    my $sl  = $seq->seq;
    $sl     =~ s/[-*]+//g; #we want ungapped lengths
    $seq->seq($sl);
    $L{$id} = $seq->length;

}

#############READ JALVIEW ANNOT
foreach my $file (@ARGV) {
    open(F, $file) or die "Cannot open JALVIEW_ANNOT_FILE! $!\n";
    while (<F>) {
        chomp; 
        my @f=split/\t/;
        
        unless (@f == 6) {
            warn "Num fields != 6! Is this a jalview annotation file?\n" ;
            next;
        }
        
        my $f = Bio::Graphics::Feature->new(-start=>$f[3], -end=>$f[4], -display_name=>$f[5]);
        my @temp            = $motifs{$f[1]} ? @{$motifs{$f[1]}} : () ;
        @temp               = sort {$a->start <=> $b->start} @temp,$f;
        @{$motifs{$f[1]}}   = @temp;
        $k{$f[5]}++;
        $dom_freq{$f[5]}++;
        push @{$ML{$f[5]}}, $f->length;
    }
}

#############DOMAIN LENGTHS
if ($opts_l) { #get average motif/domain length
    foreach my $m (keys %ML) {
        my ($c,$sum);
        foreach (@{$ML{$m}}) {
            $sum+=$_;
            $c++;
        }
        $ML{$m}=$sum/$c;
    }
}
else {
    undef %ML;
}

##############DOMAIN FREQUENCY PENALTIES
my $seq_count 		= scalar keys %L;
my $min_freq;
open(F, ">$opts_o/$opts_n.domain_penalties.txt") or die "Cannot create $opts_o/$opts_n.domain_penalties.txt! $!\n";
foreach my $d (keys (%dom_freq)) {
	$dom_freq{$d}  = $opts_l ? $dom_freq{$d}/($seq_count*$ML{$d}) : $dom_freq{$d}/$seq_count;
	$min_freq = $dom_freq{$d} if !$min_freq || $dom_freq{$d}< $min_freq;
    print F sprintf("%s\t%.4f\n", $d, $dom_freq{$d});
}
close F;



##############MAX DOMAIN LENGTH
my $max_length;
foreach my $p (keys %motifs) {
	$max_length=scalar @{$motifs{$p}} if !$max_length || @{$motifs{$p}}>$max_length ;
}

my $min_score	= ($gap_open * $max_length) / $max_length;
my $max_score	= ($max_length * (1+(1-$min_freq)))/$max_length;

##############DOMAIN ALIGN
my (@a,@b,@alignment,%seen,@matrix,@rows,@cols,$ic,$ec,@IDs,%K,%MAX,%CLUSTERS);
$ic =0;
@IDs=sort keys %L;
foreach my $i (@IDs) {
		push @cols, $i unless $seen{$i};
		$seen{$i}++;

		@a	= map{$_->display_name} @{$motifs{$i}};
		@a	= screen(@a);
		
        $ec = 0;
        foreach my $e (@IDs) {
			push @cols, $e unless $seen{$e};

			@b	= map{$_->display_name} @{$motifs{$e}};
			@b	= screen(@b);
			my ($mm)= sort {$b<=>$a} scalar(@a), scalar(@b);

			my %u;
			map {$u{$_}{a}++} @a;
			map {$u{$_}{b}++} @b;
			my $common;
			foreach (keys %u) {
				$common++ if $u{$_}{a} && $u{$_}{b};
			}
			unless ($common) { 
				# seqs don't have common domains
				my $score = $mm ? ($gap_open * $mm) / $mm : $min_score; 
                $matrix[$ic][$ec]=1-normalize_score($score);
                if ($i ne $e & (!$MAX{$i} || $MAX{$i}{score}< normalize_score($score))) {
                    $MAX{$i}{score}=normalize_score($score);
                    $MAX{$i}{seq}=$e;
                }

            }
            else {
                my $simple = Algorithm::NeedlemanWunsch->new(\&simple_scheme);
                if ($gap_open < $gap_extension) {
                    $simple->gap_open_penalty($gap_open);
                    $simple->gap_extend_penalty($gap_extension);
                }
                @alignment = ();
                my $score = $simple->align(\@a,\@b,	{
                    align => \&prepend_align,
                    shift_a => \&prepend_first_only,
                    shift_b => \&prepend_second_only ,
                    select_align => \&postpone_gap
                });

                my $cur = $score/scalar @alignment;
				$matrix[$ic][$ec]=$cur ? 1-normalize_score($cur) : 1-normalize_score($min_score);
                if ($i ne $e & (!$MAX{$i} || $MAX{$i}{score}< normalize_score($cur))) {
                    $MAX{$i}{score}=normalize_score($cur);
                    $MAX{$i}{seq}=$e;
                }
        }
        $ec++;
        $seen{$e}++;
        $seen{$e."||".$i}++;
        $seen{$i."||".$e}++;
    }
    $ic++;
}

##############EXPORT MATRIX
my $mat = Bio::Matrix::Generic->new(-colnames=>\@IDs, -rownames=>\@IDs, -values=>\@matrix);
print $P join("\t", $mat->column_names),"\n";
foreach my $n ($mat->row_names) {
	print $P join("\t", $n,$mat->row($n)),"\n";
}

##############CLUSTERING
if ($opts_k) {
	$r->run("library(fpc)");
	$r->send(qq{m<-as.matrix(read.table("$PN",sep="\t",as.is=TRUE,header=TRUE,row.names=1))});
	$r->send(qq{k<-pamkCBI(m,krange=$opts_k:1,criterion="$criterion",usepam=TRUE,scaling=FALSE,diss=TRUE)});

    my ($clusters) 	= $r->get(q{as.numeric(k$partition)});
	my ($names)		= $r->get(q{names(k$partition)});
    open(K,">$opts_o/$opts_n.clusters.txt") or die "Cannot create $opts_o/$opts_n.clusters.txt! $!\n";
    for (my $j=0; $j<@$clusters; $j++) {
        $K{$names->[$j]}=$clusters->[$j];
        print K join("\t", $names->[$j], "cluster_".$clusters->[$j]),"\n";
        push @{$CLUSTERS{$clusters->[$j]}},$names->[$j];
    }
	print "Clustered ", scalar @IDs, " sequences according to $criterion into ", $r->get(q{k$nc}), " clusters\n";
	$r->send(q{rm(m,k)});
}

##############GRAPHICAL OUTPUT
#define graphic colors by domain order
@ord = sort {$a cmp $b} keys %k;
my $i=0;
map {$k{$_}=++$i} @ord;
undef $i;

foreach my $cl (sort keys %CLUSTERS) {
    draw_profiles($cl, $CLUSTERS{$cl},($opts_p ? \*OUT : "NO"),($opts_r ? \*REPORT : "NO"));
}

exit;

###############SUBROUTINES

sub draw_profiles {
    my ($cluster,$ids,$OUT,$REPORT) = @_;

    my ($max);
    foreach (@$ids) {
        $max=$L{$_} if !$max || $max < $L{$_};
    }
    my (%SEEN,@out);
    foreach my $k (sort {$MAX{$b}{score}<=>$MAX{$a}{score} } @$ids) {
        push @out, $k unless $SEEN{$k};
        my $b=$MAX{$k}{seq};
        push @out, $b if !$SEEN{$b} && grep {$_ eq $b} @{$CLUSTERS{$cluster}} ;
        print $OUT join("\t", "cluster_".$cluster, $k, $MAX{$k}{seq}, $MAX{$k}{score}),"\n" if $OUT ne "NO" && (!$SEEN{$k} || !$SEEN{b});
        $SEEN{$k}++;
        $SEEN{$b}++;
    }
    die "ERROR\n". join(" ", @out). "\n". join(" ", @$ids)."\n" unless scalar @out == scalar @$ids;

    my $panel	= Bio::Graphics::Panel->new(-length    => $max ? $max : 1000,
        -width     		=> 800,
        -pad_left  		=> 50,
        -pad_right 		=> 10,
        -key_style 		=> 'left',
        -image_class	=>'GD::SVG',
        -spacing		=> 0,
        -empty_tracks	=> "dashed",
        -auto_pad       => 1,
    );
    my @colors=Bio::Graphics::Panel->color_names;
    while (@colors < @ord) { warn "Too many distinct domains/motifs! Colors will be duplicated!\n"; push @colors, @colors;}

    my (%mot_counts);
    foreach my $p (@out) {
        my $full = Bio::Graphics::Feature->new(-display_name => $p,-start=>1,-end=>$L{$p});
        $panel->add_track($full,
            -glyph   => 'arrow',
            -tick    => 2,
            -fgcolor => 'black',
            -double  => 1,
            -label=>1,
        );
        my $track = $panel->add_track(
            -glyph          => 'generic',
            -height		    => 7,
            -label          => 1,
            -label_spacing  => 0,
            -label_font     => 'gdTinyFont',
            -label_position	=> 'top',
        );

        foreach my $m (@{$motifs{$p}}) {
            $m->add_tag_value("bgcolor",$colors[$k{$m->display_name}]);
            $track->add_feature($m);
            $mot_counts{$m->display_name}++;
        }
    }
    if ($opts_r) {
        my $i=0;

        my @mots  = sort keys %mot_counts;
        $r->send(q{p<-NULL;});
        foreach my $m (@mots) {
            $i++;
            my $n       = $mot_counts{$m};
            my $freq    = sprintf("%.2f", $n/scalar(@out));
            my $N       = scalar(@out);
            $r->send(qq{ Q<-chisq.test(as.table(cbind(c($n,$N),c(($N*$opts_r),$N))),simulate=TRUE) });
            $r->send(qq{ p[$i]<-Q\$p.value });
        }
        $r->send(q{q<-p.adjust(p,method="BH")});
        my $P   = $r->get("p");
        my $Q   = $r->get("q");
        foreach my $p (@out) {
            $i=0;
            foreach my $m (@mots) {
                my $pmotc = scalar(grep({$m eq $_->display_name} @{$motifs{$p}}));
                my $n       = $mot_counts{$m};
                my $freq    = sprintf("%.2f", $n/scalar(@out));
                my $N       = scalar(@out);
                my $E       = sprintf("%.2f",$N*$opts_r);

                if ($pmotc && $Q->[$i]< $opts_q) { # significantly outlying
                    print $REPORT join("\t", "cluster_".$cluster,$p, $m,$pmotc, $n, $E, $N, $freq, sprintf("%.2f",$E/$N), $P->[$i],$Q->[$i], $pmotc >0 ? "PRESENT" :"ABSENT" ),"\n";
                }
                elsif (! $pmotc && $Q->[$i] >= $opts_a) { # not outlying and not found in seq
                    print $REPORT join("\t", "cluster_".$cluster,$p, $m, $pmotc,$n, $E, $N, $freq, sprintf("%.2f",$E/$N), $P->[$i],$Q->[$i], $pmotc >0 ? "PRESENT" :"ABSENT" ),"\n";
                }
                $i++;
            }
        }
    }
    open(G, ">$opts_o/$opts_n.cluster_$cluster.svg") or die $!;
    print G $panel->svg();

}

sub print_scores {
	my ($i,$e, $score, $shared,$file) = @_;
	my $norm 	= (($score- $min_score) / ($max_score - $min_score))*(1-0)+0;

    unless ( $e eq $i || $seen{$e."||".$i} || $seen{$e."||".$i}) {
        print $file join ("\t", $i, $e, $shared, sprintf("%.${comma}g",$score), sprintf("%.${comma}g",$norm) ),"\n";
    }
}
sub normalize_score {
	my ($score) = @_;
	my $norm 	= (($score- $min_score) / ($max_score - $min_score))*(1-0)+0;
	return($norm);
}


sub simple_scheme {
	if (!@_) {
		return 0;
	}
	return ($_[0] eq $_[1]) ? 1+(1-$dom_freq{$_[0]}) : -1;
}

sub prepend_align {
	my ($i, $j) = @_;
	unshift @alignment, [$a[$i], $b[$j]];
}

sub prepend_first_only {
	my $i = shift;
	unshift @alignment, [$a[$i], undef];
}

sub prepend_second_only {
	my $j = shift;
	unshift @alignment, [undef, $b[$j]];
}

sub postpone_gap {
	my $arg = shift;

	if (exists($arg->{shift_a})) {
		prepend_first_only($arg->{shift_a});
		return 'shift_a';
	} elsif (exists($arg->{shift_b})) {
		prepend_second_only($arg->{shift_b});
		return 'shift_b';
	} else {
		prepend_align(@{$arg->{align}});
		return 'align';
	}
}
sub screen {
	my (@in) = @_;
	if (@opts_i) {
		my @tt;
		foreach my $d (@in) { push @tt, $d if ! grep({$d =~ /$_/} @opts_i);}
		@in=@tt;
	}
	return @in;
}
