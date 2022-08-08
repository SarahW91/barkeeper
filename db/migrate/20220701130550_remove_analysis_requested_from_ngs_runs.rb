class RemoveAnalysisRequestedFromNgsRuns < ActiveRecord::Migration[5.2]
  def change
    remove_column :ngs_runs, :analysis_requested
  end
end
