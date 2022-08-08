class ChangeDeviseModuleFields < ActiveRecord::Migration[5.2]
  def change
    remove_column :users, :reset_password_token
    remove_column :users, :reset_password_sent_at

    add_column :users, :failed_attempts, :integer, :default => 0
    add_column :users, :locked_at, :datetime
  end
end
