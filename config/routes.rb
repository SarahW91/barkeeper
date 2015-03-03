GBOLapp::Application.routes.draw do

  resources :primer_reads

  root :to => "home#search"

  match 'help',    to: 'home#help',    via: 'get'
  match 'about',   to: 'home#about',   via: 'get'
  match 'contact', to: 'home#contact', via: 'get'
  match 'search', to: 'home#search', via: 'get'

  resources :issues

  resources :higher_order_taxons

  resources :species do
    collection do
      post :import
      get :filter
    end
    member do
      get 'show_individuals'
    end
  end

  resources :primer_pos_on_genomes

  resources :alignments

  resources :projects

  resources :news

  resources :shelves

  resources :labs

  resources :freezers

  resources :lab_racks

  resources :marker_sequences

  resources :orders

  resources :individuals do
    collection do
      get :filter
    end
  end

  resources :tissues

  resources :statuses

  resources :primers do
    collection do
      post :import
    end
  end

  resources :primer_reads do
    collection do
      post :import
      post :batch_create
    end
    member do
      get 'trim'
      get 'assign'
      get 'reverse'
      get 'restore'
      get 'fasta'
      post 'change_base'
    end
  end

  resources :plant_plates

  resources :micronic_plates

  resources :markers

  resources :isolates do
    collection do
      get 'filter'
      post :import
    end
  end

  resources :families do
    collection do
      get 'filter'
    end
    member do
      get 'show_species'
    end
  end


  resources :contigs do
    collection do
      get 'filter'
      get 'assemble_all'
      get 'pde_all'
    end
    member do
      get 'pde'
      get 'fasta'
      get 'fasta_trimmed'
      get 'fasta_raw'
      get 'overlap'
    end
  end

  #hack: avoid malicious users to directly type in the sign-up route
  #later: use authorization system to
  devise_scope :user do
    get "/users/sign_up",  :to => "home#about"
  end

  devise_for :users, :controllers => {:registrations => "registrations"}
  resources :users

  require 'sidekiq/web'
  mount Sidekiq::Web => '/sidekiq'
end