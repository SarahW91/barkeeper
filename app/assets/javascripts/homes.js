/*
 * BarKeeper - A versatile web framework to assemble, analyze and manage DNA
 * barcoding data and metadata.
 * Copyright (C) 2022 Kai Müller <kaimueller@uni-muenster.de>, Sarah Wiechers
 * <sarah.wiechers@uni-muenster.de>
 *
 * This file is part of BarKeeper.
 *
 * BarKeeper is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * BarKeeper is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with BarKeeper.  If not, see <http://www.gnu.org/licenses/>.
 */

jQuery(function() {
    $.ajax({
        type: "GET",
        contentType: "application/json; charset=utf-8",
        url: 'homes/1/background_image_urls', // This will break if there is ever more than the one "Home" record
        dataType: 'json',
        success: function (data) {
            var images = [];

            for(var i in data) {
                images.push([data[i]]);
            }

            rotateBackgrounds(data);
        },
        error: function (_result) {
            console.error("Error getting data.");
        }
    });

    function rotateBackgrounds(background_image_urls) {
        var index = 0;
        var image_div = $('#background-helper-div');

        setInterval(function () {
            image_div.animate({opacity: 0}, 500, function () {
                image_div.css('background-image', 'url(' + background_image_urls[index] + ')');
                index++;
                image_div.animate({opacity: 1}, 1000, function () {
                    if (index == background_image_urls.length) index = 0;
                });
            });
        }, 5000);
    }
});
