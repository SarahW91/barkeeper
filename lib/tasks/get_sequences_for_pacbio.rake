namespace :data do
  task get_sequences_for_pacbio: :environment do
    name_list = "GBoL2691 GBoL2694 GBoL2696 GBoL2711 GBoL2718 GBoL2695 GBoL2724 GBoL2726 GBoL2729 GBoL2731 GBoL2732
GBoL2740 GBoL2744 GBoL2748 GBoL2770 GBoL2776 GBoL2759 GBoL2787 GBoL2815 GBoL2819 GBoL2790 GBoL2846 GBoL2834 GBoL2593
GBoL2630 GBoL2631 GBoL2641 GBoL2649 GBoL2658 GBoL2666 GBoL2672 GBoL2682 GBoL3650 GBoL3652 GBoL3658 GBoL3661 GBoL3662
GBoL3697 GBoL3699 GBoL3700 GBoL3703 GBoL3704 GBoL3705 GBoL3727 GBoL3707 GBoL3701 GBoL3739 GBoL3740 GBoL3457 GBoL3467
GBoL3459 GBoL3474 GBoL3472 GBoL3476 GBoL3477 GBoL3481 GBoL3482 GBoL3505 GBoL3517 GBoL3484 GBoL3521 GBoL3522 GBoL3523
GBoL3535 GBoL3537 GBoL3554 GBoL3568 GBoL3569 GBoL3573 GBoL3571 GBoL3588 GBoL3589 GBoL3590 GBoL3591 GBoL3592 GBoL3593
GBoL3594 GBoL3595 GBoL3596 GBoL3597 GBoL3603 GBoL3613 GBoL3615 GBoL3616 GBoL3617 GBoL3618 GBoL3619 GBoL3641 GBoL3625
GBoL3626 GBoL3627 GBoL3635 GBoL3636 GBoL3606 GBoL3637 GBoL3541".split

    Marker.gbol_marker.each do |marker|
      name_list_marker = name_list.collect { |name| name + "_#{marker.name}" }
      marker_sequences = MarkerSequence.in_project(5).order(:name).where(name: name_list_marker)

      File.open("#{marker.name}_sanger.fasta", 'w') do |file|
        file.puts MarkerSequenceSearch.fasta(marker_sequences, metadata: false)
      end

      File.open("#{marker.name}_taxonomy.tax", 'w') do |file|
        file.puts MarkerSequenceSearch.taxonomy_file(marker_sequences)
      end
    end
  end
end
