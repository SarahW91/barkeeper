class ContigDatatable

  #ToDo: fig out if this inclusion is necessary. Found on https://gist.github.com/jhjguxin/4544826, but unclear if makes sense. "delegate" statement alone does not work.
  include Rails.application.routes.url_helpers

  delegate :url_helpers, to: 'Rails.application.routes'
  delegate :params, :link_to, :h, to: :@view

  # def initialize(view, need_verify, imported, duplicates)
  def initialize(view, contigs_to_show)
    @view = view
    @contigs_to_show = contigs_to_show
  end

  def as_json(options = {})
    {
        sEcho: params[:sEcho].to_i,
        iTotalRecords: Contig.count,
        iTotalDisplayRecords: contigs.total_entries,
        aaData: data
    }
  end

  private

  def data
    contigs.map do |contig|

      assembled='No'
      if contig.assembled
        assembled='Yes'
      end

      species_name=''
      species_id=0
      individual_name=''
      individual_id=0

      if contig.try(:isolate).try(:individual).try(:species)
        species_name=contig.isolate.individual.species.name_for_display
        species_id=contig.isolate.individual.species.id
        individual_name=contig.isolate.individual.specimen_id
        individual_id=contig.isolate.individual.id
      end

      [
          # check_box_tag("contig_ids[]", contig.id),
          link_to(contig.name, edit_contig_path(contig)),
          link_to(species_name, edit_species_path(species_id)),
          link_to(individual_name, edit_individual_path(individual_id)),
          assembled,
          contig.updated_at.in_time_zone("CET").strftime("%Y-%m-%d %H:%M:%S"),
          link_to('Delete', contig, method: :delete, data: { confirm: 'Are you sure?' })
      ]
    end
  end

  def contigs
    @contigs ||= fetch_contigs
  end

  def fetch_contigs

    case @contigs_to_show

      when "duplicates"
        names_with_multiple = Contig.group(:name).having("count(name) > 1").count.keys
        contigs=Contig.where(name: names_with_multiple).order("#{sort_column} #{sort_direction}")
      when "need_verification"
        contigs = Contig.need_verification.order("#{sort_column} #{sort_direction}")
      when "imported"
        contigs=Contig.externally_edited.order("#{sort_column} #{sort_direction}")
      when "caryophyllales_verified"
        contigs=Contig.caryo_matK.verified.order("#{sort_column} #{sort_direction}")
      when "caryophyllales_need_verification"
        contigs=Contig.caryo_matK.need_verification.order("#{sort_column} #{sort_direction}")
      when "caryophyllales_not_assembled"
        contigs=Contig.caryo_matK.not_assembled.order("#{sort_column} #{sort_direction}")
      when "festuca_verified"
        contigs=Contig.festuca.verified.order("#{sort_column} #{sort_direction}")
      when "festuca_need_verification"
        contigs=Contig.festuca.need_verification.order("#{sort_column} #{sort_direction}")
      when "festuca_not_assembled"
        contigs=Contig.festuca.not_assembled.order("#{sort_column} #{sort_direction}")
      else
        contigs = Contig.order("#{sort_column} #{sort_direction}")
    end

    contigs = contigs.page(page).per_page(per_page)

    if params[:sSearch].present?
      contigs = contigs.where("contigs.name ILIKE :search", search: "%#{params[:sSearch]}%")
    end
    contigs
  end

  def page
    params[:iDisplayStart].to_i/per_page + 1
  end

  def per_page
    params[:iDisplayLength].to_i > 0 ? params[:iDisplayLength].to_i : 10
  end

  def sort_column
    columns = %w[name species_id individual_id assembled updated_at]
    columns[params[:iSortCol_0].to_i]
  end

  def sort_direction
    params[:sSortDir_0] == "desc" ? "desc" : "asc"
  end
end