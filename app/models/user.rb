#
# BarKeeper - A versatile web framework to assemble, analyze and manage DNA
# barcoding data and metadata.
# Copyright (C) 2022 Kai Müller <kaimueller@uni-muenster.de>, Sarah Wiechers
# <sarah.wiechers@uni-muenster.de>
#
# This file is part of BarKeeper.
#
# BarKeeper is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# BarKeeper is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with BarKeeper.  If not, see <http://www.gnu.org/licenses/>.
#

# frozen_string_literal: true

class User < ApplicationRecord
  include ProjectRecord

  devise :database_authenticatable, :registerable, :rememberable, :timeoutable, :lockable, :validatable

  belongs_to :lab
  has_many :contig_searches, dependent: :destroy
  has_many :individual_searches, dependent: :destroy
  has_many :marker_sequence_searches, dependent: :destroy

  validates_presence_of :email # Necessary for Devise

  enum role: %i[guest user supervisor admin]
  enum responsibility: %i[taxonomy lab bulk_delete_contigs]
end
