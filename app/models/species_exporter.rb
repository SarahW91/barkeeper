#
# Barcode Workflow Manager - A web framework to assemble, analyze and manage DNA
# barcode data and metadata.
# Copyright (C) 2020 Kai Müller <kaimueller@uni-muenster.de>, Sarah Wiechers
# <sarah.wiechers@uni-muenster.de>
#
# This file is part of Barcode Workflow Manager.
#
# Barcode Workflow Manager is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# Barcode Workflow Manager is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with Barcode Workflow Manager.  If not, see
# <http://www.gnu.org/licenses/>.
#
# frozen_string_literal: true

# Write SPECIMENS & STATUS to Excel-XML (xls) for use by ZFMK for their "Portal / db : bolgermany.de "
class SpeciesExporter < ApplicationRecord
  has_one_attached :species_export
  validates :species_export, content_type: [:xml]

  def create_species_export(project_id)
    file_to_upload = File.open('species_export.xls', 'w')

    file_to_upload.write(xml_string(project_id))
    file_to_upload.close

    self.species_export.attach(io: File.open('species_export.xls'), filename: 'species_export.xls', content_type: 'application/xml')
    save!

    # Remove local file after uploading it to S3
    File.delete('species_export.xls') if File.exist?('species_export.xls')
  end

  def xml_string(project_id)
    @species = Species.includes(family: { order: :higher_order_taxon }).in_project(project_id)

    @header_cells = ['UAbteilung/Klasse',
                     'Ordnung',
                     'Familie',
                     'GBoL5_TaxID',
                     'Gattung',
                     'Art',
                     'Autor',
                     'Subspecies/Varietät',
                     'Autor (ssp./var.)']

    builder = Nokogiri::XML::Builder.new do |xml|
      xml.comment("Generated by gbol5.de web app on #{Time.zone.now}")
      xml.Workbook('xmlns' => 'urn:schemas-microsoft-com:office:spreadsheet',
                   'xmlns:o' => 'urn:schemas-microsoft-com:office:office',
                   'xmlns:x' => 'urn:schemas-microsoft-com:office:excel',
                   'xmlns:ss' => 'urn:schemas-microsoft-com:office:spreadsheet',
                   'xmlns:html' => 'https://www.w3.org/TR/REC-html40') do
        xml.Worksheet('ss:Name' => 'Sheet1') do
          xml.Table do
            # Header
            xml.Row do
              @header_cells.each do |o|
                xml.Cell do
                  xml.Data('ss:Type' => 'String') do
                    xml.text(o)
                  end
                end
              end
            end

            @species.each do |species|
              xml.Row do
                # Subdivision/Class
                xml.Cell do
                  xml.Data('ss:Type' => 'String') do
                    xml.text(species.try(:family).try(:order).try(:higher_order_taxon).try(:name))
                  end
                end

                # Order
                xml.Cell do
                  xml.Data('ss:Type' => 'String') do
                    xml.text(species.try(:family).try(:order).try(:name))
                  end
                end

                # Family
                xml.Cell do
                  xml.Data('ss:Type' => 'String') do
                    xml.text(species.try(:family).try(:name))
                  end
                end

                # GBoL5_TaxID
                xml.Cell do
                  xml.Data('ss:Type' => 'String') do
                    xml.text(species.id)
                  end
                end

                # Genus
                xml.Cell do
                  xml.Data('ss:Type' => 'String') do
                    xml.text(species.genus_name)
                  end
                end

                # Species
                xml.Cell do
                  xml.Data('ss:Type' => 'String') do
                    xml.text(species.species_epithet)
                  end
                end

                # Author
                xml.Cell do
                  xml.Data('ss:Type' => 'String') do
                    xml.text(species.author)
                  end
                end

                # Subspecies/Variety
                xml.Cell do
                  xml.Data('ss:Type' => 'String') do
                    xml.text(species.infraspecific)
                  end
                end

                # Author (ssp./var.)
                xml.Cell do
                  xml.Data('ss:Type' => 'String') do
                    xml.text(species.author_infra)
                  end
                end
              end
            end
          end
        end
      end
    end

    builder.to_xml
  end
end
