# frozen_string_literal: true

class AddAlignedConsToContigs < ActiveRecord::Migration
  def change
    add_column :contigs, :aligned_cons, :string
  end
end
