# frozen_string_literal: true

class IndividualsController < ApplicationController
  include ProjectConcern

  load_and_authorize_resource

  before_action :set_individual, only: %i[show edit update destroy]

  def index
    respond_to do |format|
      format.html
      format.json { render json: IndividualDatatable.new(view_context, nil, current_project_id) }
    end
  end

  def export_as_csv
    authorize! :export_as_csv, :individual

    send_data(Individual.to_csv(current_project_id),
              filename: "specimen_project_#{Project.find(current_project_id).name}.csv",
              type: 'application/csv')
  end

  def problematic_specimens
    # Liste aller Bundesländer to check 'state_province' against:
    @states = %w[Baden-Württemberg Bayern Berlin Brandenburg Bremen Hamburg Hessen Mecklenburg-Vorpommern Niedersachsen Nordrhein-Westfalen Rheinland-Pfalz Saarland Sachsen Sachsen-Anhalt Schleswig-Holstein Thüringen]
    @individuals = []

    Individual.in_project(current_project_id).each do |i|
      if i.country == 'Germany'
        @individuals.push(i) unless @states.include? i.state_province
      end
    end
  end

  def filter
    @individuals = Individual.where('individuals.specimen_id ilike ?', "%#{params[:term]}%").in_project(current_project_id)
    render json: @individuals.map{ |individual| {:id=> individual.id, :name => individual.specimen_id }}
  end

  def show; end

  def new
    @individual = Individual.new
  end

  def edit; end

  def create
    @individual = Individual.new(individual_params)
    @individual.add_project(current_project_id)

    respond_to do |format|
      if @individual.save
        format.html { redirect_to individuals_path, notice: 'Individual was successfully created.' }
        format.json { render :show, status: :created, location: @individual }
      else
        format.html { render :new }
        format.json { render json: @individual.errors, status: :unprocessable_entity }
      end
    end
  end

  def update
    respond_to do |format|
      if @individual.update(individual_params)
        format.html { redirect_to individuals_path, notice: 'Individual was successfully updated.' }
        format.json { render :show, status: :ok, location: @individual }
      else
        format.html { render :edit }
        format.json { render json: @individual.errors, status: :unprocessable_entity }
      end
    end
  end

  def destroy
    @individual.destroy
    respond_to do |format|
      format.html { redirect_to individuals_url }
      format.json { head :no_content }
    end
  end

  private

  # Use callbacks to share common setup or constraints between actions.
  def set_individual
    @individual = Individual.find(params[:id])
  end

  # Never trust parameters from the scary internet, only allow the white list through.
  def individual_params
    params.require(:individual).permit(:specimen_id, :DNA_bank_id, :collector, :specimen_id, :herbarium_code,
                                       :herbarium_id, :country, :state_province, :locality, :latitude, :longitude,
                                       :latitude_original, :longitude_original, :elevation, :exposition, :habitat,
                                       :substrate, :life_form, :collectors_field_number, :collected, :determination,
                                       :revision, :confirmation, :comments, :taxon_id, :taxon_name, :tissue_id,
                                       project_ids: [])
  end
end
