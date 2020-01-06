# frozen_string_literal: true

# Write SPECIMENS & STATUS to Excel-XML (xls) for use by ZFMK for their "Portal / db : bolgermany.de "
class SpecimenExporter < ApplicationRecord
  include ActionView::Helpers

  has_attached_file :specimen_export,
                    path: '/specimens_export.xls'

  # Validate content type
  validates_attachment_content_type :specimen_export, content_type: /\Aapplication\/xml/
  # Validate filename
  validates_attachment_file_name :specimen_export, matches: [/xls\Z/]

  def create_specimen_export(project_id)
    file_to_upload = File.open('specimens_export.xls', 'w')
    xml_file(file_to_upload, project_id)
    file_to_upload.close

    self.specimen_export = File.open('specimens_export.xls')
    save!
  end

  def xml_file(file, project_id)
    markers = Marker.gbol_marker

    @states = %w[Baden-Württemberg Bayern Berlin Brandenburg Bremen Hamburg
                 Hessen Mecklenburg-Vorpommern Niedersachsen Nordrhein-Westfalen
                 Rheinland-Pfalz Saarland Sachsen Sachsen-Anhalt Schleswig-Holstein
                 Thüringen]

    @header_cells = ['GBOL5 specimen ID', 'Feldnummer', 'Institut', 'Sammlungs-Nr.',
                     'Familie', 'Taxon-Name', 'Erstbeschreiber Jahr', 'evtl. Bemerkung Taxonomie',
                     'Name', 'Datum', 'Gewebetyp und Menge', 'Anzahl Individuen', 'Fixierungsmethode',
                     'Entwicklungsstadium', 'Sex', 'evtl. Bemerkungen zur Probe', 'Fundortbeschreibung',
                     'Region', 'Bundesland', 'Land', 'Datum', 'Sammelmethode', 'Breitengrad',
                     'Längengrad', 'Benutzte Methode', 'Ungenauigkeitsangabe', 'Höhe/Tiefe [m]',
                     'Habitat', 'Sammler', 'Feld-Sammelnummer', 'Behörde']

    markers.each do |m|
      @header_cells << "#{m.name} - URL"
      @header_cells << "#{m.name} - Marker Sequence"
      @header_cells << "#{m.name} - Genbank ID"
    end

    xml = ::Builder::XmlMarkup.new(target: file, indent: 2)

    xml.instruct! :xml, version: '1.0', encoding: 'UTF-8'
    xml.comment!("Generated by gbol5.de web app on #{Time.zone.now}")

    xml.Workbook('xmlns' => 'urn:schemas-microsoft-com:office:spreadsheet',
                 'xmlns:o' => 'urn:schemas-microsoft-com:office:office',
                 'xmlns:x' => 'urn:schemas-microsoft-com:office:excel',
                 'xmlns:ss' => 'urn:schemas-microsoft-com:office:spreadsheet',
                 'xmlns:html' => 'https://www.w3.org/TR/REC-html40') do
      xml.Worksheet('ss:Name' => 'Sheet1') do
        xml.Table do
          xml.Row do
            @header_cells.each { |header| cell_tag(xml, header) }
          end

          # Individuals in current project
          Individual.includes(species: :family).in_project(project_id).find_each do |individual|
            xml.Row do
              cell_tag(xml, individual.id) # GBOL5 specimen ID

              # Feldnummer
              if individual.collectors_field_number&.include?('s.n.') || individual.collectors_field_number&.include?('s. n.')
                cell_tag(xml, '')
              else
                cell_tag(xml, individual.collectors_field_number)
              end

              cell_tag(xml, individual.herbarium_code) # Institut

              # Sammlungsnummer
              if individual.specimen_id == '<no info available in DNA Bank>'
                cell_tag(xml, '')
              else
                cell_tag(xml, individual.specimen_id)
              end

              cell_tag(xml, individual.try(:species).try(:family).try(:name)) # Familie
              cell_tag(xml, individual.try(:species).try(:name_for_display)) # Taxonname
              cell_tag(xml, individual.try(:species).try(:author)) # Erstbeschreiber Jahr
              cell_tag(xml, '') # evtl. Bemerkung Taxonomie
              cell_tag(xml, individual.determination) # Name
              cell_tag(xml, '') # Datum
              cell_tag(xml, 'Blattmaterial') # Gewebetyp und Menge
              cell_tag(xml, '') # Anzahl Individuen
              cell_tag(xml, 'Silica gel') # Fixierungsmethode
              cell_tag(xml, '') # Entwicklungsstadium
              cell_tag(xml, '') # Sex
              cell_tag(xml, "gbol5.de/individuals/#{individual.id}/edit") # evtl. Bemerkungen zur Probe
              cell_tag(xml, individual.locality) # Fundortbeschreibung
              cell_tag(xml, '') # Region

              # Bundesland
              if individual.country == 'Germany' || individual.country == 'Deutschland'
                # Tests first if is a Bundesland; outputs nothing if other crap was entered in this field:
                if @states.include? individual.state_province
                  cell_tag(xml, individual.state_province)
                else
                  cell_tag(xml, '')
                end
              else
                cell_tag(xml, 'Europa') # Stuff from Schweiz etc.
              end

              cell_tag(xml, individual.country) # Land
              cell_tag(xml, individual.collection_date) # Datum
              cell_tag(xml, '') # Sammelmethode
              cell_tag(xml, number_with_precision(individual.latitude, precision: 5)) # Breitengrad
              cell_tag(xml, number_with_precision(individual.longitude, precision: 5)) # Längengrad
              cell_tag(xml, '') # Benutzte Methode
              cell_tag(xml, '') # Ungenauigkeitsangabe
              cell_tag(xml, individual.elevation) # Höhe/Tiefe [m]
              cell_tag(xml, individual.habitat) # Habitat
              cell_tag(xml, individual.collector) # Sammler
              cell_tag(xml, individual.collectors_field_number) # Feld-Sammelnummer
              cell_tag(xml, '') # Behörde

              markers.each do |marker|
                # Find longest marker sequence per GBoL marker
                longest_sequence = MarkerSequence.joins(isolate: :individual).includes(:contigs)
                                                 .where('individuals.id = ?', individual.id)
                                                 .where('marker_id = ?', marker.id).limit(1)
                                                 .group(:id).order('MAX(CHAR_LENGTH(sequence)) desc').first

                # URL zum Contig in GBOL5 WebApp
                url = longest_sequence&.contigs&.first&.id ? "gbol5.de/contigs/#{longest_sequence&.contigs&.first&.id}/edit" : ''
                cell_tag(xml, url)

                cell_tag(xml, longest_sequence&.sequence) # Markersequenz
                cell_tag(xml, longest_sequence&.genbank) # Genbank ID
              end
            end
          end
        end
      end
    end
  end

  private

  def cell_tag(xml, string)
    string ||= ''
    xml.Cell do
      xml.Data(string, 'ss:Type' => 'String')
    end
  end
end
