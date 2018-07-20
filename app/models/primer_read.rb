class PrimerRead < ApplicationRecord
  include ProjectRecord

  belongs_to :contig
  belongs_to :partial_con
  belongs_to :primer
  has_many :issues

  has_attached_file :chromatogram,
                    :default_url => '/chromatograms/primer_read.scf'

  # Do_not_validate_attachment_file_type :chromatogram

  # Validate content type
  validates_attachment_content_type :chromatogram, :content_type => /\Aapplication\/octet-stream/

  # Validate filename
  validates_attachment_file_name :chromatogram, :matches => [/scf\Z/, /ab1\Z/]

  before_create :default_name

  scope :assembled, -> {use_for_assembly.where(:assembled => true)}
  scope :not_assembled, -> {use_for_assembly.where(:assembled => false)}

  scope :use_for_assembly, ->  { trimmed.where(:used_for_con => true)}
  scope :not_used_for_assembly, ->  { trimmed.where(:used_for_con => false)}

  scope :trimmed, -> { where.not(:trimmedReadStart => nil)}
  scope :not_trimmed, -> { where(:trimmedReadStart => nil)}

  scope :processed, -> {where(:processed => true)}
  scope :unprocessed, -> {where(:processed => false)}
  scope :contig_not_verified, -> {joins(:contig).where(:contigs => {:verified => false, :verified_by => nil})}

  validates_attachment_presence :chromatogram

  def self.in_higher_order_taxon(higher_order_taxon_id)
    count = 0

    HigherOrderTaxon.find(higher_order_taxon_id).orders.each do |ord|
      ord.families.each do |fam|
        fam.species.each do  |sp|
          sp.individuals.each do |ind|
            ind.isolates.each do |iso|
              iso.contigs.each do |con|
                count += con.primer_reads.count
              end
            end
          end
        end
      end
    end

    count
  end

  def file_name_id
    self.name.gsub('.', "_#{self.id}.")
  end

  def slice_to_json(start_pos, end_pos)
    # get trace position corresponding to first / last aligned_peaks index that exists and is not -1:

    start_pos_trace=start_pos
    end_pos_trace=end_pos

    xstart=self.aligned_peak_indices[start_pos_trace]

    while xstart == -1
      start_pos_trace+=1
      xstart=self.aligned_peak_indices[start_pos_trace]
    end

    # does aligned_peaks index  exist?
    if self.aligned_peak_indices[end_pos_trace]
      xend=self.aligned_peak_indices[end_pos_trace]
      #else use last:
    else
      end_pos_trace=self.aligned_peak_indices.count-1
      xend=self.aligned_peak_indices[end_pos_trace]
    end

    # find first that isnt -1:
    while xend == -1 and end_pos_trace > 0
      end_pos_trace-=1
      xend=self.aligned_peak_indices[end_pos_trace]
    end

    # create hash with x-pos as key, ya, yc, … as value
    traces=Hash.new

    if xstart and xend #account for situations where nothing from this read is seen in respective contig slice/page:

      xstart-=10
      xend+=10

      (xstart..xend).each do |x|
        traces[x] = {
            :ay => self.atrace[x],
            :cy => self.ctrace[x],
            :gy => self.gtrace[x],
            :ty => self.ttrace[x]
        }
      end
    end

    # get position in original, non-trimmed, non-aligned primer read:


    # return json

    {
        :id => self.id.as_json,
        :name => self.name.as_json,
        :aligned_seq => self.aligned_seq[start_pos..end_pos].as_json,
        :aligned_qualities => self.aligned_qualities[start_pos..end_pos].as_json,
        :traces => traces.as_json,
        :aligned_peak_indices => self.aligned_peak_indices[start_pos..end_pos].as_json,
        :trimmedReadStart => self.trimmedReadStart.as_json,
        :trimmedReadEnd => self.trimmedReadEnd.as_json,
        :original_positions => self.original_positions[start_pos..end_pos].as_json
    }
  end

  def original_positions
    original_positions=Array.new

    i=self.trimmedReadStart

    self.aligned_qualities.each do |aq|

      if aq == -1
        original_positions << -1
      else
        original_positions << i
        i+=1
      end

    end

    original_positions
  end

  def contig_name
    contig.try(:name)
  end

  def contig_name=(name)
    if name == ''
      self.contig = nil
    else
      self.contig = Contig.find_or_create_by(:name => name) if name.present?
    end
  end

  def seq_for_display
    "#{self.sequence[0..30]...self.sequence[-30..-1]}" if self.sequence.present?
  end

  def trimmed_seq_for_display
    "#{self.sequence[self.trimmedReadStart..self.trimmedReadStart+30]...self.sequence[self.trimmedReadEnd-30..self.trimmedReadEnd]}" if self.sequence.present?
  end

  def name_for_display
    if self.name.length > 25
      "#{self.name[0..11]...self.name[-11..-5]}"
    else
      "#{self.name[0..-5]}"
    end
  end

  def default_name
    self.name ||= self.chromatogram.original_filename
  end

  def auto_assign
    output_message = nil
    create_issue = false

    # Try to find matching primer
    regex_read_name = /^([A-Za-z0-9]+)(.*)_([A-Za-z0-9-]+)\.(scf|ab1)$/ # match group 1: GBoL number, 2: stuff, 3: primer name, 4: file extension
    name_components = self.name.match(regex_read_name)

    if name_components
      primer_name = name_components[3]
      name_variants_t7 = %w[T7promoter T7 T7-1] # T7 is always forward
      name_variants_m13 = %w[M13R-pUC M13-RP M13-RP-1] #M13R-pU

      #logic if T7promoter or M13R-pUC.scf:
      if name_variants_t7.include? primer_name
        rgx = /^_([A-Za-z0-9]+)_([A-Za-z0-9]+)$/
        matches = name_components[2].match(rgx) # --> uv2
        primer_name = matches[1]

        primer = Primer.where("primers.name ILIKE ?", "#{primer_name}").first
        primer ||= Primer.where("primers.alt_name ILIKE ?", "#{primer_name}").first

        if primer
          primer_name = matches[2] if primer.reverse
        else
          output_message = "Cannot find primer with name #{primer_name}."
          create_issue = true
        end
      elsif name_variants_m13.include? primer_name
        rgx = /^_([A-Za-z0-9]+)_([A-Za-z0-9]+)$/
        matches = name_components[2].match(rgx) # --> 4
        primer_name = matches[2] # --> uv4

        primer = Primer.where("primers.name ILIKE ?", "#{primer_name}").first
        primer ||= Primer.where("primers.alt_name ILIKE ?", "#{primer_name}").first

        if primer
          primer_name = matches[1] unless primer.reverse
        else
          output_message = "Cannot find primer with name #{primer_name}."
          create_issue = true
        end
      else
        # Leave primer_name as is
      end

      # Find & assign primer

      primer = Primer.where("primers.name ILIKE ?", "#{primer_name}").first
      primer ||= Primer.where("primers.alt_name ILIKE ?", "#{primer_name}").first

      if primer
        self.update(:primer_id => primer.id, :reverse => primer.reverse)

        # find marker that primer belongs to
        marker = primer.marker

        if marker
          # Try to find matching isolate
          isolate_component = name_components[1] # GBoL number

          # BGBM cases:
          regex_db_number = /^.*(DB)[\s_]?([0-9]+)(.*)_([A-Za-z0-9-]+)\.(scf|ab1)$/ # match group 1: DNABank number, 2: stuff, 3: primer name, 4: file extension
          db_number_name_components = self.name.match(regex_db_number)

          if db_number_name_components
            isolate_component = "#{db_number_name_components[1]} #{db_number_name_components[2]}" # DNABank number
          end

          isolate = Isolate.where("isolates.lab_nr ILIKE ?", isolate_component.to_s).first
          isolate ||= Isolate.where("isolates.dna_bank_id ILIKE ?", isolate_component.to_s).first
          isolate ||= Isolate.create(:lab_nr => isolate_component)

          if db_number_name_components
            isolate.update(:dna_bank_id => isolate_component)
          end

          self.update(:isolate_id => isolate.id)

          # Figure out which contig to assign to
          matching_contig = Contig.where("contigs.marker_id = ? AND contigs.isolate_id = ?", marker.id, isolate.id).first

          if matching_contig
            self.contig = matching_contig
            self.save
            output_message = "Assigned to contig #{matching_contig.name}."
          else
            # Create new contig, auto assign to primer, copy, auto-name
            contig = Contig.new(:marker_id => marker.id, :isolate_id => isolate.id, :assembled => false)

            contig.generate_name
            contig.save

            self.contig = contig
            self.save

            output_message = "Created contig #{contig.name} and assigned primer read to it."
          end
        else
          output_message = "No marker assigned to primer #{primer.name}."
          create_issue = true
        end
      else
        output_message = "Cannot find primer with name #{primer_name}."
        create_issue = true
      end
    else
      output_message = "No match for #{regex_read_name} in name."
      create_issue = true
    end

    if create_issue
      i = Issue.create(:title => output_message, :primer_read_id => self.id)
    else # Everything worked
      self.contig.update(:assembled => false, :assembly_tried => false)
    end

    { :msg => output_message, :create_issue => create_issue }
  end

  def get_position_in_marker(p)
    # get position in marker

    pp=nil

    if self.trimmed_seq
      if self.reverse
        pp= p.position-self.trimmed_seq.length
      else
        pp= p.position
      end
    end

    pp
  end

  def copy_to_db

  end

  def auto_trim(write_to_db)

    msg=nil
    create_issue = false

    #get local copy from s3

    dest = Tempfile.new(self.chromatogram_file_name)
    dest.binmode
    self.chromatogram.copy_to_local_file(:original, dest.path)

    begin

      chromatogram_ff1 = nil
      p = /\.ab1$/

      if self.chromatogram_file_name.match(p)
        chromatogram_ff1 = Bio::Abif.open(dest.path)
      else
        chromatogram_ff1 = Bio::Scf.open(dest.path)
      end

      chromatogram1 = chromatogram_ff1.next_entry

      if self.reverse
        chromatogram1.complement!()
      end

      sequence = chromatogram1.sequence.upcase

      #copy chromatogram over into db
      if self.sequence.nil? or write_to_db
        self.update(:sequence => sequence)
        self.update(:base_count => sequence.length)
      end
      if self.qualities.nil? or write_to_db
        self.update(:qualities => chromatogram1.qualities)
      end
      if self.atrace.nil? or write_to_db
        self.update(:atrace => chromatogram1.atrace)
        self.update(:ctrace => chromatogram1.ctrace)
        self.update(:gtrace => chromatogram1.gtrace)
        self.update(:ttrace => chromatogram1.ttrace)
        self.update(:peak_indices => chromatogram1.peak_indices)
      end

      se = Array.new

      se = self.trim_seq(chromatogram1.qualities, self.min_quality_score, self.window_size, self.count_in_window)

      #se = self.trim_seq_inverse(chromatogram1.qualities)

      if se
        if se[0]>=se[1] # trimming has not found any stretch of bases > min_score
          msg='Quality too low - no stretch of readable bases found.'
          create_issue=true
          self.update(:used_for_con => false)
        elsif se[2] > 0.6

          msg="Quality too low - #{(se[2]*100).round}% low-quality base calls in trimmed sequence."
          create_issue=true

          self.update(:used_for_con => false)
        else

          # everything works:

          self.update(:trimmedReadStart => se[0]+1, :trimmedReadEnd => se[1]+1, :used_for_con => true)
          #:position => self.get_position_in_marker(self.primer)
          msg='Sequence trimmed.'
        end
      else
        msg='Quality too low - no stretch of readable bases found.'
        create_issue=true
        self.update(:used_for_con => false)
      end
    rescue
      msg='Sequence could not be trimmed - no scf/ab1 file or no quality scores?'
      create_issue = true
      self.update(:used_for_con => false)
    end

    if create_issue
      i=Issue.create(:title => msg, :primer_read_id => self.id)
      self.update(:used_for_con=>false)
    end

    {:msg => msg, :create_issue => create_issue}

  end

  def trimmed_seq
    if trimmedReadStart.nil? or trimmedReadEnd.nil?
      nil
    else
      if trimmedReadEnd > trimmedReadStart
        self.sequence[(self.trimmedReadStart-1)..(self.trimmedReadEnd-1)] if self.sequence.present?
        #cleaned_sequence = raw_sequence.gsub('-', '') # in case basecalls in pherogram have already '-' - as in some crappy seq. I got from BN
      else
        nil
      end
    end
  end

  def trimmed_quals
    if trimmedReadStart.nil? or trimmedReadEnd.nil?
      nil
    else
      if trimmedReadEnd > trimmedReadStart
        self.qualities[(self.trimmedReadStart-1)..(self.trimmedReadEnd-1)] if self.qualities.present?
      else
        nil
      end
    end
  end

  def trimmed_and_cleaned_seq
    self.trimmed_seq.upcase.gsub /[^ACTGN-]+/, 'N'
  end

  def get_aligned_peak_indices

    if self.trimmedReadStart

      aligned_peak_indices = Array.new

      pi=self.trimmedReadStart-2

      if self.aligned_qualities

        self.aligned_qualities.each do |aq|
          if aq==-1
            aligned_peak_indices << -1
          else
            aligned_peak_indices << self.peak_indices[pi]
            pi+=1
          end
        end

        self.update_columns(aligned_peak_indices: aligned_peak_indices)

      end
    end

  end


  #deactivated cause even worse than original
  # def trim_seq_inverse(qualities)
  #
  #   #settings
  #   min_quality_score = 20
  #   c=16
  #   t=20
  #
  #   #final coordinates:
  #
  #   trimmed_read_start = 0
  #   trimmed_read_end = qualities.length
  #
  #   #intermediate coordinates:
  #
  #   trimmed_read_start1 = 0
  #   trimmed_read_end1 = qualities.length
  #
  #   trimmed_read_start2 = qualities.length
  #   trimmed_read_end2 = 0
  #
  #   # --- find readstart:
  #
  #   for i in 0..qualities.length-t
  #     #extract window of size t
  #
  #     count=0
  #
  #     for k in i...i+t
  #       if qualities[k]>=min_quality_score
  #         count += 1
  #       end
  #     end
  #
  #     if count>=c
  #       trimmed_read_start1 = i
  #       break
  #     end
  #
  #   end
  #
  #   # stop when already at seq end
  #   if i >= qualities.length-t
  #     return nil
  #   end
  #
  #   # -- find read-end1, BUT THIS TIME COMING FROM SAME DIRECTION:
  #   #asking: When does quality stop to obey the above quality condition?
  #
  #   for i in trimmed_read_start1..qualities.length-t
  #     #extract window of size t
  #
  #     count=0
  #
  #     for k in i...i+t
  #       if qualities[k]>=min_quality_score
  #         count += 1
  #       end
  #     end
  #
  #     if count<c
  #       trimmed_read_end1 = i
  #       break
  #     end
  #
  #   end
  #
  #   #### same from other dir:
  #
  #   # --- find readend2:
  #
  #   i =qualities.length
  #
  #   while i > 0
  #     #extract window of size t
  #
  #     # k=i
  #     count=0
  #
  #     for k in i-t...i
  #       if qualities[k]>=min_quality_score
  #         count += 1
  #       end
  #     end
  #
  #     if count>=c
  #       trimmed_read_end2 = i
  #       break
  #     end
  #
  #     i-=1
  #
  #   end
  #
  #   # -- find read_start2, BUT THIS TIME COMING FROM SAME DIRECTION:
  #   #asking: When does quality stop to obey the above quality condition?
  #
  #   i =trimmed_read_end2
  #
  #   while i > 0
  #     #extract window of size t
  #
  #     # k=i
  #     count=0
  #
  #     for k in i-t...i
  #       if qualities[k]>=min_quality_score
  #         count += 1
  #       end
  #     end
  #
  #     if count<c
  #       trimmed_read_start2 = i
  #       break
  #     end
  #
  #     i-=1
  #
  #   end
  #
  #
  #   # choose longer fragment
  #
  #   puts trimmed_read_start1
  #   puts trimmed_read_end1
  #
  #   puts trimmed_read_start2
  #   puts trimmed_read_end2
  #
  #
  #
  #   if (trimmed_read_end2-trimmed_read_start2) > (trimmed_read_end1-trimmed_read_start1)
  #     trimmed_read_end=trimmed_read_end2
  #     trimmed_read_start=trimmed_read_start2
  #   else
  #     trimmed_read_end=trimmed_read_end1
  #     trimmed_read_start=trimmed_read_start1
  #   end
  #
  #   #check if xy% < min_score:
  #   ctr_bad=0
  #   ctr_total=0
  #   for j in trimmed_read_start...trimmed_read_end
  #     if qualities[j]<min_quality_score
  #       ctr_bad+=1
  #     end
  #     ctr_total+=1
  #   end
  #
  #   [trimmed_read_start, trimmed_read_end, ctr_bad.to_f/ctr_total.to_f]
  #
  #
  #
  #
  # end


  # old version

  def trim_seq(qualities, min_quality_score, t, c)

    trimmed_read_start = 0
    trimmed_read_end = qualities.length

    # --- find readstart:

    for i in 0..qualities.length-t
      #extract window of size t

      count=0

      for k in i...i+t
        if qualities[k]>=min_quality_score
          count += 1
        end
      end

      if count>=c
        trimmed_read_start = i
        break
      end

    end

    #fine-tune:  if bad bases are at beginning of last window, cut further until current base's score >= min_qual:

    ctr=trimmed_read_start

    for a in ctr..ctr+t
      if qualities[a]>=min_quality_score
        trimmed_read_start = a
        break
      end
    end


    # --- find readend:

    i =qualities.length
    while i > 0
      #extract window of size t

      # k=i
      count=0

      for k in i-t...i
        if qualities[k]>=min_quality_score
          count += 1
        end
      end

      if count>=c
        trimmed_read_end = i
        break
      end

      i-=1

    end

    #fine-tune:  if bad bases are at beginning of last window, go back until current base's score >= min_qual:

    while i > trimmed_read_end-t
      if qualities[i]>=min_quality_score
        break
      end
      i-=1
    end
    trimmed_read_end = i

    #check if xy% < min_score:
    ctr_bad=0
    ctr_total=0
    for j in trimmed_read_start...trimmed_read_end
      if qualities[j]<min_quality_score
        ctr_bad+=1
      end
      ctr_total+=1
    end

    [trimmed_read_start, trimmed_read_end, ctr_bad.to_f/ctr_total.to_f]

  end


end