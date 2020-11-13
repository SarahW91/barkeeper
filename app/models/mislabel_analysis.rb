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

class MislabelAnalysis < ApplicationRecord
  belongs_to :marker
  has_many :mislabels, dependent: :destroy
  has_one :marker_sequence_search, dependent: :destroy
  has_and_belongs_to_many :marker_sequences

  def percentage_of_mislabels
    ((mislabels.size / total_seq_number.to_f) * 100).round(2) if total_seq_number&.positive?
  end

  def analyse_on_server
    Net::SSH.start('xylocalyx.uni-muenster.de', 'kai', keys: ['/home/sarah/.ssh/gbol_xylocalyx']) do |session|
      # Check if SATIVA.sh is already running
      running = session.exec!("pgrep -f \"SATIVA.sh\"")

      if running.empty?
        sequences = "#{Rails.root}/#{title}.fasta"
        tax_file = "#{Rails.root}/#{title}.tax"

        analysis_dir = "/data/data1/sarah/SATIVA/#{title}"

        File.open(sequences, 'w+') do |f|
          f.write(marker_sequence_search.analysis_fasta(true))
        end

        File.open(tax_file, 'w+') do |f|
          f.write(marker_sequence_search.taxon_file(true))
        end

        # Create analysis directory
        session.exec!("mkdir #{analysis_dir}")

        # Upload analysis input files
        session.scp.upload! tax_file, analysis_dir
        session.scp.upload! sequences, analysis_dir

        # Delete local files
        FileUtils.rm(sequences)
        FileUtils.rm(tax_file)

        # Start analysis on server
        session.exec!("/home/kai/analysis-scripts/SATIVA.sh #{title} #{id}")
      end
    end
  end

  # Check if recent SATIVA results exist and download them
  def download_results
    exists = false
    results_remote = "/data/data1/sarah/SATIVA/#{title}/#{title}.mis"
    results = "#{Rails.root}/#{title}.mis"

    # Check if file exists before download
    Net::SFTP.start('xylocalyx.uni-muenster.de', 'kai', keys: ['/home/sarah/.ssh/gbol_xylocalyx']) do |sftp|
      sftp.stat(results_remote) do |response|
        if response.ok?
          # Download result file
          sftp.download!(results_remote, results)
          exists = true
        end
      end
    end

    # Import analysis result and then remove them
    if exists
      import(results)
      FileUtils.rm(results)

      # Also remove last automated analysis for this marker from web app
      last_analysis = MislabelAnalysis.where(automatic: true, marker: marker).order(created_at: :desc).last
      last_analysis.destroy unless last_analysis.id == id # Avoid self-destruct!

      return true
    else
      return false
    end
  end

  private

  # Import SATIVA result table (*.mis)
  def import(file)
    column_names = %w[SeqID MislabeledLevel OriginalLabel ProposedLabel
                      Confidence OriginalTaxonomyPath ProposedTaxonomyPath
                      PerRankConfidence]

    File.open(file, 'r').each do |row|
      next if row.start_with?(';')

      row = Hash[[column_names, row.split("\t")].transpose]

      mislabel = Mislabel.new
      mislabel.level = row['MislabeledLevel']
      mislabel.confidence = row['Confidence'].to_f
      mislabel.proposed_label = row['ProposedLabel']
      mislabel.proposed_path = row['ProposedTaxonomyPath']
      mislabel.path_confidence = row['PerRankConfidence']
      mislabel.save

      name = row['SeqID'].split('_')[0..1].join('_')
      marker_sequence = MarkerSequence.find_by_name(name)
      if marker_sequence.nil? && name.start_with?('DB')
        marker_sequence = MarkerSequence.find_by_name(name.gsub('DB', 'DB '))
      end

      if marker_sequence
        marker_sequence.mislabels << mislabel

        mislabels << mislabel
        marker_sequences << marker_sequence
      end
    end
  end
end
