class FamiliesController < ApplicationController
  include ProjectConcern

  load_and_authorize_resource

  before_action :set_family, only: [:show, :edit, :update, :destroy]

  # GET /families
  # GET /families.json

  def index
    @families = Family.includes(:order).in_project(current_project_id).order('name asc')
    respond_to :html, :json
  end

  def filter
    @families = Family.in_project(current_project_id).order(:name).where("families.name ilike ?", "%#{params[:term]}%")
    render json: @families.map(&:name)
  end

  def show_species
    respond_to do |format|
      format.html
      format.json { render json: SpeciesDatatable.new(view_context, params[:id], nil, current_project_id) }
    end
  end

  # GET /families/1
  # GET /families/1.json
  def show
  end

  # GET /families/new
  def new
    @family = Family.new
  end

  # GET /families/1/edit
  def edit
  end

  # POST /families
  # POST /families.json
  def create
    @family = Family.new(family_params)
    @family.add_project(current_project_id)

    respond_to do |format|
      if @family.save
        format.html { redirect_to families_path, notice: 'Family was successfully created.' }
        format.json { render :show, status: :created, location: @family }
      else
        format.html { render :new }
        format.json { render json: @family.errors, status: :unprocessable_entity }
      end
    end
  end

  # PATCH/PUT /families/1
  # PATCH/PUT /families/1.json
  def update
    respond_to do |format|
      if @family.update(family_params)
        format.html { redirect_to families_path, notice: 'Family was successfully updated.' }
        format.json { render :show, status: :ok, location: @family }
      else
        format.html { render :edit }
        format.json { render json: @family.errors, status: :unprocessable_entity }
      end
    end
  end

  # DELETE /families/1
  # DELETE /families/1.json
  def destroy
    @family.destroy
    respond_to do |format|
      format.html { redirect_to families_url }
      format.json { head :no_content }
    end
  end

  private

  # Use callbacks to share common setup or constraints between actions.
  def set_family
    @family = Family.find(params[:id])
  end

  # Never trust parameters from the scary internet, only allow the white list through.
  def family_params
    params.require(:family).permit(:name, :author, :order_id, :term, :project_ids => [])
  end
end
