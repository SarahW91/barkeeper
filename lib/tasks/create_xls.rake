namespace :data do

  desc 'export specimen info to ZFMK readable xls format'
  task :create_xls => :environment do
    SpecimenExport.perform_async(current_project_id)
  end

end