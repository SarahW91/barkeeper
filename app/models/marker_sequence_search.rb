# frozen_string_literal: true

class MarkerSequenceSearch < ApplicationRecord
  include Export

  belongs_to :user
  belongs_to :project
  belongs_to :mislabel_analysis

  validates :title, uniqueness: { :scope=> :user_id }

  enum has_warnings: %i[both yes no]

  def marker_sequences
    @marker_sequences ||= find_marker_sequences
  end

  def analysis_fasta(include_singletons, metadata=false, label_warnings=false)
    sequences = include_singletons ? marker_sequences : remove_singletons(marker_sequences)

    MarkerSequenceSearch.fasta(sequences, { warnings: label_warnings, metadata: metadata })
  end

  def taxon_file(include_singletons)
    sequences = include_singletons ? marker_sequences : remove_singletons(marker_sequences)
    MarkerSequenceSearch.taxonomy_file(sequences)
  end

  private

  def find_marker_sequences
    if project_id
      marker_sequences = MarkerSequence.in_project(project_id).order(:name)
    else
      marker_sequences = MarkerSequence.order(:name)
    end

    if name.present?
      if name.include?(',')
        names = name.split(',')
        marker_sequences = marker_sequences.where(name: names)
      else
        marker_sequences = marker_sequences.where('marker_sequences.name ilike ?', "%#{name}%")
      end
    end

    if verified != 'both'
      marker_sequences = marker_sequences.verified if verified == 'verified'
      marker_sequences = marker_sequences.not_verified if verified == 'unverified'
    end

    marker_sequences = marker_sequences.has_species if has_species.present?

    marker_sequences = marker_sequences.no_isolate if no_isolate.present?

    if has_warnings != 'both'
      marker_sequences = marker_sequences.unsolved_warnings if has_warnings == 'yes'
      marker_sequences = marker_sequences.where.not(id: MarkerSequence.unsolved_warnings.pluck(:id)) if has_warnings == 'no'
    end

    marker_sequences = marker_sequences.joins(:marker).where('markers.name ilike ?', "%#{marker}%") if marker.present?

    marker_sequences = marker_sequences.joins(isolate: { individual: { species: { family: { order: :higher_order_taxon } } } }).where('higher_order_taxa.name ilike ?', "%#{higher_order_taxon}%") if higher_order_taxon.present?

    marker_sequences = marker_sequences.joins(isolate: { individual: { species: { family: :order } } }).where('orders.name ilike ?', "%#{self.order}%") if self.order.present?

    marker_sequences = marker_sequences.joins(isolate: { individual: { species: :family } }).where('families.name ilike ?', "%#{family}%") if family.present?

    marker_sequences = marker_sequences.joins(isolate: { individual: :species }).where('species.composed_name ilike ?', "%#{species}%") if species.present?

    marker_sequences = marker_sequences.joins(isolate: :individual).where('individuals.specimen_id ilike ?', "%#{specimen}%") if specimen.present?

    marker_sequences = marker_sequences.joins(:contigs).where('contigs.verified_by = ?', User.find_by_name(verified_by)&.id) if verified_by.present?

    marker_sequences = marker_sequences.where('length(marker_sequences.sequence) >= ?', min_length) if min_length.present?
    marker_sequences = marker_sequences.where('length(marker_sequences.sequence) <= ?', max_length) if max_length.present?

    marker_sequences = marker_sequences.where('marker_sequences.created_at >= ?', min_age.midnight) if min_age.present?
    marker_sequences = marker_sequences.where('marker_sequences.created_at <= ?', max_age.end_of_day) if max_age.present?

    marker_sequences = marker_sequences.where('marker_sequences.updated_at >= ?', min_update.midnight) if min_update.present?
    marker_sequences = marker_sequences.where('marker_sequences.updated_at <= ?', max_update.end_of_day) if max_update.present?

    marker_sequences
  end

  def remove_singletons(sequences)
    cnt = sequences.joins(isolate: [individual: :species]).distinct.reorder('species.species_component').group('species.species_component').count
    sequences.where(species: { species_component: cnt.select { |_, v| v > 1 }.keys })
  end
end
