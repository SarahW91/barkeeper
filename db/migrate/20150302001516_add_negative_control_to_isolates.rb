class AddNegativeControlToIsolates < ActiveRecord::Migration
  def change
    add_column :isolates, :negative_control, :boolean, :default => false
  end
end