module Bio::Graphics::Helper
  def self.plotgene(my_panel, gene_id)
    require 'ensembl'
    Ensembl::Core::DBConnection.connect('homo_sapiens',58)
    transcripts_track = my_panel.add_track('Protein Coding Transcripts', :glyph => {'utr' => :spliced, 'cds' => :directed_spliced}, :colour => {'utr' => [1,0,0],'cds' => [0,0,1]})
    gene =  Ensembl::Core::Gene.find(gene_id, :include => {:transcripts => :exon_transcripts})
    transcripts = gene.transcripts
    transcripts = transcripts.delete_if{|x| x.biotype != 'protein_coding'}
    transcripts.each do |transcript|
      utr5 = ""
      utr3 = ""
      cds = ""
      transcript.exons.sort{|x| x.seq_region_start}.each do |exon|
        case
        when exon.seq_region_end <=  transcript.coding_region_genomic_start  && exon.seq_region_start >= transcript.seq_region_start
          utr5 += "#{exon.seq_region_start}..#{exon.seq_region_end},"
        when exon.seq_region_start <= transcript.coding_region_genomic_start && exon.seq_region_end >= transcript.coding_region_genomic_start
          utr5 += "#{exon.seq_region_start}..#{transcript.coding_region_genomic_start},"
          cds += "#{transcript.coding_region_genomic_start}..#{exon.seq_region_end},"
        when exon.seq_region_start >= transcript.coding_region_genomic_start && exon.seq_region_end <= transcript.coding_region_genomic_end
          cds += "#{exon.seq_region_start}..#{exon.seq_region_end},"
        when exon.seq_region_start <= transcript.coding_region_genomic_end && exon.seq_region_end >= transcript.coding_region_genomic_end
          cds += "#{exon.seq_region_start}..#{transcript.coding_region_genomic_end},"
          utr3 += "#{transcript.coding_region_genomic_end}..#{exon.seq_region_end},"
        when  exon.seq_region_start >= transcript.coding_region_genomic_end && exon.seq_region_end <= transcript.seq_region_end
          utr3 += "#{exon.seq_region_start}..#{exon.seq_region_end},"
        end

      end
      cds.chomp!(',')
      utr5.chomp!(',')
      utr3.chomp!(',')
      if transcript.seq_region_strand == 1
        join = "join(#{utr5 + ',' unless utr5 == ""}#{cds},#{utr3 unless utr3 == ""})"
        cdsobject = Bio::Feature.new('cds', "join(#{cds})")
        utr5object = Bio::Feature.new('utr', "join(#{utr5})")
        utr3object = Bio::Feature.new('utr', "join(#{utr3})")
      else

        join = "complement(join(#{utr3 + ',' unless utr3 == ""}#{cds},#{utr5 unless utr5 == ""}))"
        cdsobject = Bio::Feature.new('cds', "complement(join(#{cds}))")
        utr5object = Bio::Feature.new('utr', "complement(join(#{utr5}))")
        utr3object = Bio::Feature.new('utr', "complement(join(#{utr3}))")
      end

      subfeatures = []
      subfeatures << cdsobject unless cds == ""
      subfeatures << utr5object unless utr5 == ""
      subfeatures << utr3object unless utr3 == ""
      transcript2 = Bio::Feature.new(transcript.stable_id,join,[],nil,subfeatures)
      transcripts_track.add_feature(transcript2,:label => transcript.stable_id, :link => "http://www.ensembl.org/Homo_sapiens/Transcript/Summary?db=core;t=#{transcript.stable_id}")
    end
    return my_panel
  end


end