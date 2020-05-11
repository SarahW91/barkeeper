# frozen_string_literal: true

source 'https://rubygems.org'

ruby '2.4.0'

gem 'rails', '5.0.0'

# Postgres
gem 'pg', '~> 0.20.0' # Rails does not work with pg 1.0.0 (fixed in Rails 5.1.5), pg 0.21 gives constant deprecation warning: https://github.com/rails/rails/issues/29521
gem 'pg_search'

# Web server and background processing
gem 'puma', '~> 3.0'
gem 'redis', '~> 3.0' # Use Redis adapter to run Action Cable in production
gem 'sidekiq'
gem 'sidekiq-client-cli'
gem 'sidekiq-limit_fetch'
gem 'sinatra', require: false # Needed to monitor sidekiq jobs

# Asset pipeline
gem 'coffee-rails', '~> 4.2' # Use CoffeeScript for .coffee assets and views
gem 'jquery-rails'
gem 'jquery-ui-rails'
gem 'sass-rails', '~> 5.0' # Use SCSS for stylesheets
gem 'uglifier', '>= 1.3.0' # Use Uglifier as compressor for JavaScript assets

# Authentication & authorization
gem 'cancancan'
gem 'devise'

# External file storage
gem 'aws-sdk'
gem 'paperclip'

# Javascript
gem 'bootstrap-multiselect_rails' # multi select boxes using bootstrap
gem 'chosen-rails' # Javascript select boxes
gem 'jquery-datatables-rails', github: 'rweng/jquery-datatables-rails'
gem 'nested_form_fields'
gem 'jquery-fileupload-rails'
gem 'jquery-turbolinks'
gem 'select2-rails' # Integrate Select2 Javascript library

gem 'ancestry' # Self-related models
gem 'bcrypt', '~> 3.1.7', platforms: :ruby
gem 'bio' # BioRuby
gem 'bootstrap-sass'
gem 'bootstrap_progressbar'
gem 'builder'
gem 'jbuilder', '~> 2.5' # Build JSON APIs with ease. Read more: https://github.com/rails/jbuilder
gem 'mime-types'
gem 'net-scp'
gem 'net-sftp'
gem 'net-ssh'
gem 'roo-xls' # Handle excel files
gem 'rubyzip' # Handle zip files
gem 'simple_form'
gem 'slim' # TODO: Used anywhere?
gem 'sprockets-rails', require: 'sprockets/railtie'
# gem 'turbolinks' # Turbolinks makes navigating your web application faster TODO: Does not work unless first changes to js code, in particular data-tables (see http://guides.rubyonrails.org/working_with_javascript_in_rails.html#turbolinks)
gem 'whenever', require: false # Runs scheduled jobs via cron
gem 'will_paginate', '> 3.0'
gem 'will_paginate-bootstrap'

# Windows does not include zoneinfo files, so bundle the tzinfo-data gem
gem 'tzinfo-data', platforms: %i[mingw mswin x64_mingw jruby]

group :development do
  gem 'capistrano',         require: false
  gem 'capistrano-bundler', require: false
  gem 'capistrano-rails',   require: false
  gem 'capistrano-rbenv',   require: false
  gem 'capistrano-sidekiq', require: false, github: 'seuros/capistrano-sidekiq'
  gem 'capistrano3-puma', '~> 3.1', require: false
  gem 'better_errors' # Better error page for Rack apps
  gem 'bullet' # Checks for n+1 queries
  gem 'binding_of_caller' # Extends features of better_errors
  gem 'meta_request' # Supporting gem for Google Chrome Rails Panel
  # gem 'spring' # Spring speeds up development by keeping your application running in the background
  # gem 'spring-watcher-listen', '~> 2.0.0'
  gem 'yard' # Documentation generation tool
end

group :test do
  gem 'capybara', '~> 2'
  gem 'capybara-selenium'
  gem 'database_cleaner'
  gem 'factory_bot_rails'
  gem 'faker'
  gem 'guard-rspec'
  gem 'launchy'
  gem 'webdrivers'
  gem 'rails_best_practices'
  gem 'rspec-rails'
  gem 'shoulda-matchers'
  gem 'shoulda-callback-matchers'
  gem 'simplecov'
end
