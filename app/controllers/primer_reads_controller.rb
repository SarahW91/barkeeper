class PrimerReadsController < ApplicationController
  before_action :authenticate_user!, :except => [:edit, :index]

  before_action :set_primer_read, only: [:go_to_pos, :do_not_use_for_assembly, :use_for_assembly, :change_left_clip, :change_right_clip, :edit, :fasta, :reverse, :restore, :assign, :trim, :show, :update, :change_base, :destroy]

  def do_not_use_for_assembly
    @primer_read.update(:used_for_con => false, :assembled => false)
    if @primer_read.contig
      if @primer_read.contig.primer_reads.where(:used_for_con => true).count <= 4
        @primer_read.contig.auto_overlap
        msg='Assembly finished.'
      else
        ContigAssembly.perform_async(@primer_read.contig.id)
        msg='Assembly started in background.'
      end
    end
    redirect_back(fallback_location: primer_reads_path, notice: msg)
  end

  def use_for_assembly
    @primer_read.update(:used_for_con => true)
    if @primer_read.contig
      if @primer_read.contig.primer_reads.where(:used_for_con => true).count <= 4
        @primer_read.contig.auto_overlap
        msg='Assembly finished.'
      else
        ContigAssembly.perform_async(@primer_read.contig.id)
        msg='Assembly started in background.'
      end
    end
    redirect_back(fallback_location: primer_reads_path, notice: msg)
  end

  def duplicates
    respond_to do |format|
      format.html
      format.json { render json: PrimerReadDatatable.new(view_context, true)}
    end
  end

  def import
    PrimerRead.import(params[:file])
    redirect_to root_url, notice: 'Chromatogram imported.'
  end

  def reads_without_contigs
    respond_to do |format|
      format.html
      format.json { render json: PrimerReadWithoutContigDatatable.new(view_context)}
    end
  end

  # GET /primer_reads
  # GET /primer_reads.json
  def index
    respond_to do |format|
      format.html
      format.json { render json: PrimerReadDatatable.new(view_context, false)}
    end
  end

  # GET /primer_reads/1
  # GET /primer_reads/1.json
  def show
  end

  # GET /primer_reads/new
  def new
    @primer_read = PrimerRead.new
  end

=begin
  def batch_create
    if params[:chromatograms]
      params[:chromatograms].each { |chromatogram|
        PrimerRead.create(chromatogram: chromatogram, name: chromatogram.original_filename)
      }
    end
    @primer_reads = PrimerRead.all
    #render :index
    redirect_to root_url, notice: 'Chromatograms imported.'
  end
=end

  def trim

    msg_hash = @primer_read.auto_trim(false)

    if msg_hash[:create_issue]
      redirect_to edit_primer_read_path, alert: msg_hash[:msg]
    else
      redirect_to edit_primer_read_path, notice: msg_hash[:msg]
    end

  end

  # tries to extract associated primer and isolate from primer read name (in turn based on uploaded scf file name):
  def assign

    msg_hash = @primer_read.auto_assign

    if msg_hash[:create_issue]
      redirect_to edit_primer_read_path, alert: msg_hash[:msg]
    else
      redirect_to edit_primer_read_path, notice: msg_hash[:msg]
    end

  end

  def reverse
    unless @primer_read.reverse
      begin
        @primer_read.update(:reverse => true)
        @primer_read.auto_trim(true)
        redirect_to edit_primer_read_path, notice: 'Reversed.'
      rescue
        redirect_to edit_primer_read_path, alert: 'Could not reverse'
      end
    else
      redirect_to edit_primer_read_path, alert: 'Already reversed.'
    end
  end

  def restore
    if @primer_read.reverse

      begin
        @primer_read.update(:reverse => false)
        @primer_read.auto_trim(true)

        redirect_to edit_primer_read_path, notice: 'Restored non-reversed state.'
      rescue
        redirect_to edit_primer_read_path, alert: 'Could not restore.'
      end
    else
      redirect_to edit_primer_read_path, alert: 'Read is not reversed.'
    end
  end

  # GET /primer_reads/1/edit
  def edit
    #@primer_read = PrimerRead.includes(:primer, :contig).find(params[:id]) #TODO Add select statement here to initially NOT load chromatogram?
  end

  def go_to_pos
    @pos = params[:pos]
  end

  # POST /primer_reads
  # POST /primer_reads.json

  def create

    @primer_read = PrimerRead.new(primer_read_params)

    if @primer_read.save
      PherogramProcessing.perform_async(@primer_read.id)
      if @primer_read.sequence.nil?
        @primer_read.update(:sequence => '')
      end
    end

  end


  # PATCH/PUT /primer_reads/1
  # PATCH/PUT /primer_reads/1.json
  def update
    respond_to do |format|
      if @primer_read.update(primer_read_params)
        format.html { redirect_back(fallback_location: edit_primer_read_path(@primer_read), notice: 'Primer read was successfully updated.') }
        format.json { render :show, status: :ok, location: @primer_read }
      else
        format.html { render :edit }
        format.json { render json: @primer_read.errors, status: :unprocessable_entity }
      end
    end
  end

  def change_base

    sequence= @primer_read.sequence
    pos=params[:position].to_i
    base=params[:base]

    # if insertions needed: handle inserting new elements in qualities etc. arrays
    if base.length > 1

      insertions_needed=base.length-1

      qualities=@primer_read.qualities
      insertions_needed.times do
        # insert placeholder element into qualities
        qualities.insert(pos+1,-10) #-1 already taken by aligned_qualities to indicate "gap" to be drawn in contig view
      end

      peak_indices=@primer_read.peak_indices

      #compute new indices based on existing neighbors:
      left_index=peak_indices[pos]
      right_index=peak_indices[pos+1]
      distance=right_index-left_index
      x_increment=(distance/base.length)

      x=left_index
      insertions_needed.times do
        x=x+x_increment
        # insert placeholder element into peak_indices:
        peak_indices.insert(pos+1,x)
      end

    end

    #in all cases (replacement, insertion & deletion) insert string:
    sequence[pos] = base
    @primer_read.update(:sequence => sequence)

    # do not auto assemble after base change!
    # if @primer_read.contig
    #   ContigAssembly.perform_async(@primer_read.contig.id)
    # end

    head :ok
  end

  def change_left_clip
    puts @primer_read.name
    @primer_read.update(:trimmedReadStart => params[:position].to_i)

    # do not auto assemble after clipping changed!
    # if @primer_read.contig
    #   ContigAssembly.perform_async(@primer_read.contig.id)
    # end

    head :ok
  end

  def change_right_clip
    @primer_read.update(:trimmedReadEnd => params[:position].to_i)

    # do not auto assemble after clipping changed!
    # if @primer_read.contig
    #   ContigAssembly.perform_async(@primer_read.contig.id)
    # end

    head :ok
  end

  # DELETE /primer_reads/1
  # DELETE /primer_reads/1.json
  def destroy
    @primer_read.destroy
    respond_to do |format|
      format.html { redirect_back(fallback_location: primer_reads_url) }
      format.json { head :no_content }
    end
  end


  def fasta
    str=">#{@primer_read.name} \n#{@primer_read.trimmed_seq}"
    # send_data str
    send_data(str, :filename => "#{@primer_read.name}.fas", :type => "application/txt")
  end

  private
  # Use callbacks to share common setup or constraints between actions.
  def set_primer_read
    @primer_read = PrimerRead.includes(:contig => { :isolate => :individual }).find(params[:id])
  end

  # Never trust parameters from the scary internet, only allow the white list through.
  def primer_read_params
    params.require(:primer_read).permit(:pos, :overwritten, :assembled, :comment,
                                        :quality_string,
                                        :used_for_con, :file, :name, :sequence, :pherogram_url, :chromatogram, :primer_id, :contig_id,
                                        :contig_name, :isolate_id, :chromatograms, :trimmedReadEnd, :trimmedReadStart,
                                        :min_quality_score, :count_in_window, :window_size,
                                        :position, :base)
  end

end