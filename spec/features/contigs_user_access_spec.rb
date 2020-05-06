require "rails_helper"

RSpec.feature "User access to contigs", type: :feature, js: true do
  before(:each) { @general_project = Project.find_or_create_by(name: 'All') }

  scenario "Visitor can access contigs index and filter records without login" do
    visit_index_and_filter
  end

  scenario "Visitor can access contigs edit page without login" do
    can_edit
  end

  scenario "Visitor cannot update a contig" do
    cannot_update_contig
  end

  scenario "Guest user cannot update a contig" do
    user = create(:user, role: 'guest')
    sign_in user

    cannot_update_contig
  end

  scenario "User can update a contig" do
    user = create(:user, role: 'user')
    sign_in user

    can_update_contig
  end

  scenario "Admin can update a contig" do
    user = create(:user, role: 'admin')
    sign_in user

    can_update_contig
  end

  scenario "Visitor cannot delete a contig without login" do
    cannot_delete_contig
  end

  scenario "Guest user cannot delete a contig" do
    user = create(:user, role: 'guest')
    sign_in user

    cannot_delete_contig
  end

  scenario "User can delete a contig" do
    user = create(:user, role: 'user')
    sign_in user

    can_delete_contig
  end

  scenario "Supervisor can delete a contig" do
    user = create(:user, role: 'supervisor')
    sign_in user

    can_delete_contig
  end

  scenario "Admin can delete a contig" do
    user = create(:user, role: 'admin')
    sign_in user

    can_delete_contig
  end

  def visit_index_and_filter
    first_contig = FactoryBot.create(:contig_with_taxonomy, name: "first_contig")
    second_contig = FactoryBot.create(:contig_with_taxonomy, name: "second_contig")

    species = first_contig.isolate.individual.species
    individual = first_contig.isolate.individual

    visit contigs_path

    within "h3" do
      expect(page).to have_content "Contigs"
    end

    expect(page).to have_content first_contig.name
    expect(page).to have_content second_contig.name

    # Filter by contig name partial
    find(:xpath, '//input[@type="search"]').set('second')

    expect(page).to_not have_content first_contig.name
    expect(page).to have_content second_contig.name

    # Filter by species name
    find(:xpath, '//input[@type="search"]').set(species.name_for_display)

    expect(page).to have_content first_contig.name
    expect(page).to_not have_content second_contig.name

    # Filter by specimen identifier
    find(:xpath, '//input[@type="search"]').set(individual.DNA_bank_id)

    expect(page).to have_content first_contig.name
    expect(page).to_not have_content second_contig.name
  end

  def can_edit
    species = FactoryBot.create(:species_with_individuals)
    FactoryBot.create(:individual_with_isolates, species: species)
    test_contig = FactoryBot.create(:contig, name: "test_contig", isolate: Isolate.first)

    FactoryBot.create(:contig, name: "second_test_contig")

    visit contigs_path

    within "h3" do
      expect(page).to have_content "Contigs"
    end

    expect(page).to have_content "test_contig"
    expect(page).to have_content "second_test_contig"

    click_link "test_contig"

    expect(current_path).to eq edit_contig_path(test_contig)

    within "h3" do
      expect(page).to have_content "Contig test_contig"
      expect(page).to have_link(species.name_for_display, :href => edit_species_path(species))
    end
  end

  def can_delete_contig
    test_contig = FactoryBot.create(:contig, name: "first_contig")
    FactoryBot.create(:contig, name: "second_contig")

    visit contigs_path

    find(:xpath, "//a[@href='/contigs/#{test_contig.id}'][@data-method='delete']").click
    alert = page.driver.browser.switch_to.alert
    expect(alert.text).to match "Are you sure?"
    alert.accept

    expect(page).to have_content "Contigs"
    expect(page).to have_content "second_contig"
    expect(page).to_not have_content "first_contig"

    Capybara.reset_sessions!
    DatabaseCleaner.clean
  end

  def cannot_delete_contig
    FactoryBot.create(:contig, name: "test_contig")
    FactoryBot.create(:contig, name: "second_test_contig")

    visit contigs_path

    expect {
      find(:xpath, "//a[@href='/contigs/1']").click
      alert = page.driver.browser.switch_to.alert
      expect(alert.text).to match "Are you sure?"
      alert.accept
    }.to_not change(Contig,:count)

    expect(page).to have_content "Contigs"
    expect(page).to have_content "test_contig"

    DatabaseCleaner.clean
  end

  def cannot_update_contig
    species = FactoryBot.create(:species_with_individuals)
    FactoryBot.create(:individual_with_isolates, species: species)
    test_contig = FactoryBot.create(:contig, name: "test_contig", isolate: Isolate.first)

    FactoryBot.create(:contig, name: "second_test_contig")

    visit contigs_path

    click_link 'test_contig'
    click_link 'Description'

    fill_in('Name', :with => 'Changed_Test_Name')

    expect {
      click_button 'Update'
      expect(page).to have_content "You are not authorized to access this page or perform this action."
    }.to_not change(test_contig,:name)

    within "h3" do
      expect(page).to have_content "Contig test_contig"
      expect(page).to have_link(species.name_for_display, :href => edit_species_path(species))
    end

    DatabaseCleaner.clean
  end

  def can_update_contig
    species = FactoryBot.create(:species_with_individuals)
    FactoryBot.create(:individual_with_isolates, species: species)
    test_contig = FactoryBot.create(:contig, name: "test_contig", isolate: Isolate.first)

    FactoryBot.create(:contig, name: "second_test_contig")

    visit contigs_path

    click_link 'test_contig'
    click_link 'Description'

    fill_in('Name', :with => 'Changed_Test_Name')

    click_button 'Update'
    expect(page).to have_content "Contig was successfully updated."

    within "h3" do
      expect(page).to have_content "Contig Changed_Test_Name"
      expect(page).to have_link(species.name_for_display, :href => edit_species_path(species))
    end

    DatabaseCleaner.clean
  end
end