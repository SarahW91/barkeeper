class NgsRunsController < ApplicationController
  include ProjectConcern

  load_and_authorize_resource

  http_basic_authenticate_with name: Rails.application.credentials.api_user, password: Rails.application.credentials.api_password, only: [:import, :revised_tpm]
  skip_before_action :verify_authenticity_token, only: [:import, :revised_tpm]

  before_action :set_ngs_run, only: [:show, :edit, :update, :destroy, :import, :analysis_results]

  def index
    respond_to do |format|
      format.html
      format.json { render json: NgsRunDatatable.new(view_context, current_project_id) }
    end
  end

  def show
  end

  def new
    @ngs_run = NgsRun.new
  end

  def edit
  end

  def create
    @ngs_run = NgsRun.new(ngs_run_params)
    @ngs_run.add_project(current_project_id)

    unless params[:ngs_run][:tag_primer_map].blank?
      params[:ngs_run][:tag_primer_map].each do |tpm|
        @ngs_run.tag_primer_maps.build(tag_primer_map: tpm)
      end
    end

    respond_to do |format|
      if @ngs_run.save
        format.html { redirect_to ngs_runs_path, notice: 'NGS Run was successfully created.' }
        format.json { render :edit, status: :created, location: @ngs_run }
      else
        format.html { render :new }
        format.json { render json: @ngs_run.errors, status: :unprocessable_entity }
      end
    end
  end

  def update
    params[:ngs_run].delete(:set_tag_map) if params[:ngs_run][:set_tag_map].blank?

    if params[:ngs_run][:tag_primer_map].blank?
      params[:ngs_run].delete(:tag_primer_map)
    else
      params[:ngs_run][:tag_primer_map].each do |tpm|
        # Only add TPM if Package map is available or none was added before
        if !params[:ngs_run][:set_tag_map].blank? || @ngs_run.set_tag_map.attached? || @ngs_run.tag_primer_maps.size.zero?
          map = TagPrimerMap.create(tag_primer_map: tpm)
          @ngs_run.tag_primer_maps << map
        else
          map = TagPrimerMap.create(tag_primer_map: tpm)
          @ngs_run.tag_primer_maps.delete_all
          @ngs_run.tag_primer_maps << map
        end
      end
    end

    @ngs_run.remove_tag_primer_maps(params[:delete_tag_primer_maps]) unless params[:delete_tag_primer_maps].blank?

    respond_to do |format|
      if @ngs_run.update(ngs_run_params)
        format.html { redirect_to edit_ngs_run_path(@ngs_run), notice: 'NGS Run was successfully updated.' }
        format.json { render :index, status: :ok }
      else
        format.html { render :edit }
        format.json { render json: @ngs_run.errors, status: :unprocessable_entity }
      end
    end
  end

  def destroy
    @ngs_run.destroy
    respond_to do |format|
      format.html { redirect_to ngs_runs_url, notice: 'NGS Run was successfully destroyed.' }
      format.json { head :no_content }
    end
  end

  def  submit_analysis_request
    @ngs_run.submit_request(current_user, edit_ngs_run_path(@ngs_run))
    redirect_to ngs_runs_path, notice: 'Request for analysis successfully submitted. An admin is informed. Please check in later to see results.'
  end

  def start_analysis
    if @ngs_run.check_tag_primer_map_count
      isolates_not_in_db = @ngs_run.samples_exist
      if isolates_not_in_db.blank?
        running = @ngs_run.check_server_status # Check if barcoding pipe is already running

        if running.empty?
          AnalyseNGSData.perform_async(@ngs_run.id)
          redirect_to ngs_runs_path, notice: 'Analysing data. This may take a while. Check in later to see results.'
        else
          redirect_to ngs_runs_path, alert: "Analysis could not be started: Server is busy. Try again later."
        end
      else
        redirect_back(fallback_location: ngs_runs_path,
                      alert: "Please create database entries for the following sample IDs before starting the import: #{isolates_not_in_db.join(', ')}")
      end
    else
      redirect_back(fallback_location: ngs_runs_path, alert: 'Please make sure you uploaded as many properly formatted tag primer maps as stated in the packages fasta.')
    end
  end

  def import
    render plain: "Results stored in #{params[:results_path]} will be imported in the background.\n"
    NGSResultsImporter.perform_async(@ngs_run.id, params[:results_path])
  end

  def analysis_results
    respond_to do |format|
      format.html
      format.json { render json: NgsRunResultDatatable.new(view_context, params[:id]) }
    end
  end

  def revised_tpm
    tpm = TagPrimerMap.create(tag_primer_map: params[:tpm])

    if tpm.check_tag_primer_map
      send_data(tpm.revised_tag_primer_map([5]), filename: tpm.tag_primer_map.filename, type: 'application/txt')
    else
      render plain: "The Tag Primer Map is not properly formatted. Please try again after checking the file!\n"
    end
  end

  private
  # Use callbacks to share common setup or constraints between actions.
  def set_ngs_run
    @ngs_run = NgsRun.find(params[:id])
  end

  # Never trust parameters from the scary internet, only allow the white list through.
  def ngs_run_params
    params.require(:ngs_run).permit(:name, :primer_mismatches, :quality_threshold, :tag_mismatches, :fastq_location,
                                    :set_tag_map, :higher_order_taxon_id, :delete_set_tag_map, :results)
  end
end
