class AddChildrenCountToTaxa < ActiveRecord::Migration[5.2]
  def change
    add_column :taxa, :children_count, :integer
  end
end
