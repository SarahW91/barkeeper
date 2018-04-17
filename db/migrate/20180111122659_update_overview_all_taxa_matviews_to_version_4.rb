class UpdateOverviewAllTaxaMatviewsToVersion4 < ActiveRecord::Migration[5.0]
  def change
    update_view :overview_all_taxa_matviews,
                version: 4,
                revert_to_version: 3,
                materialized: true
  end
end