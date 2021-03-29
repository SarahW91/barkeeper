# frozen_string_literal: true

class IndividualDatatable
  # TODO: fig out if this inclusion is necessary. Found on https://gist.github.com/jhjguxin/4544826, but unclear if makes sense. "delegate" statement alone does not work.

  include Rails.application.routes.url_helpers
  delegate :url_helpers, to: 'Rails.application.routes'

  delegate :params, :link_to, :h, to: :@view

  def initialize(view, taxon_id, current_default_project)
    @view = view
    @taxon_id = taxon_id
    @current_default_project = current_default_project
  end

  def as_json(_options = {})
    {
      sEcho: params[:sEcho].to_i,
      iTotalRecords: Individual.in_project(@current_default_project).count,
      iTotalDisplayRecords: individuals.total_entries,
      aaData: data
    }
  end

  private

  def data
    individuals.map do |individual|
      taxon = ''

      taxon = link_to individual.taxon.scientific_name, edit_taxon_path(individual.taxon) if individual.taxon
      herbarium = link_to individual.herbarium.acronym, edit_herbarium_path(individual.herbarium) if individual.herbarium

      [
        link_to(individual.specimen_id, edit_individual_path(individual)),
        taxon,
        herbarium,
        individual.collector,
        individual.collectors_field_number,
        individual.updated_at.in_time_zone('CET').strftime('%Y-%m-%d %H:%M:%S'),
        link_to('Delete', individual, method: :delete, data: { confirm: 'Are you sure?' })
      ]
    end
  end

  def individuals
    @individuals ||= fetch_individuals
  end

  def fetch_individuals
    if @taxon_id
      individuals = Individual.includes(:taxon, :herbarium).where(taxon_id: @taxon_id).in_project(@current_default_project).order("#{sort_column} #{sort_direction}")
    else
      individuals = Individual.includes(:taxon, :herbarium).in_project(@current_default_project).order("#{sort_column} #{sort_direction}")
    end

    individuals = individuals.page(page).per_page(per_page)

    if params[:sSearch].present?
      individuals = individuals.where('individuals.specimen_id ILIKE :search
OR taxa.scientific_name ILIKE :search
OR herbaria.acronym ILIKE :search
OR individuals.collector ILIKE :search
OR individuals.collectors_field_number ILIKE :search', search: "%#{params[:sSearch]}%")
      .references(:taxon)
    end

    individuals
  end

  def page
    params[:iDisplayStart].to_i / per_page + 1
  end

  def per_page
    params[:iDisplayLength].to_i.positive? ? params[:iDisplayLength].to_i : 10
  end

  def sort_column
    columns = %w[individuals.specimen_id taxa.scientific_name herbaria.acronym individuals.collector individuals.collectors_field_number individuals.updated_at]
    columns[params[:iSortCol_0].to_i]
  end

  def sort_direction
    params[:sSortDir_0] == 'desc' ? 'desc' : 'asc'
  end
end
