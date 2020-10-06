class NgsRunDatatable

  include Rails.application.routes.url_helpers
  delegate :url_helpers, to: 'Rails.application.routes'

  delegate :params, :link_to, :h, to: :@view


  def initialize(view, current_default_project)
    @view = view
    @current_default_project = current_default_project
  end

  def as_json(options = {})
    {
        sEcho: params[:sEcho].to_i,
        iTotalRecords: NgsRun.count,
        iTotalDisplayRecords: ngs_runs.total_entries,
        aaData: data
    }
  end

  private

  def data
    ngs_runs.map do |ngs_run|

      ngs_run_path = if ngs_run.results.attached?
                       link_to(ngs_run.name, ngs_run_path(ngs_run, anchor: "results"))
                     else
                       link_to(ngs_run.name, edit_ngs_run_path(ngs_run))
                     end

      [
        ngs_run_path,
        ngs_run.updated_at.in_time_zone("CET").strftime("%Y-%m-%d %H:%M:%S"),
        link_to('Delete', ngs_run, method: :delete, data: { confirm: 'Are you sure?' })
      ]

    end

  end

  def ngs_runs
    @ngs_runs ||= fetch_ngs_runs
  end

  def fetch_ngs_runs
    ngs_runs = NgsRun.in_project(@current_default_project).order("#{sort_column} #{sort_direction}")

    ngs_runs = ngs_runs.page(page).per_page(per_page)

    if params[:sSearch].present?
      ngs_runs = ngs_runs.where("ngs_runs.name ILIKE :search", search: "%#{params[:sSearch]}%")
    end

    ngs_runs
  end

  def page
    params[:iDisplayStart].to_i/per_page + 1
  end

  def per_page
    params[:iDisplayLength].to_i > 0 ? params[:iDisplayLength].to_i : 10
  end

  def sort_column
    columns = %w[ngs_runs.name updated_at]
    columns[params[:iSortCol_0].to_i]
  end

  def sort_direction
    params[:sSortDir_0] == "desc" ? "desc" : "asc"
  end

end