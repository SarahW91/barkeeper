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
    var body = $('#page_body');
    var helper_div = $("#background-helper-div");

    if (helper_div != null) {
        var home_id = helper_div.data('home-id');

        $.ajax({
            type: "GET",
            contentType: "application/json; charset=utf-8",
            url: 'homes/' + home_id + '/background_image_urls',
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
    }
    else {
        body.css('background-color', '#ededed');
    }

    function rotateBackgrounds(background_image_urls) {
        var index = 0;
        var image_div = $('#background-helper-div');

        body.css('background-color', '#222');

        setInterval(function () {
            image_div.animate({opacity: 0}, 600, function () {
                image_div.css('background-image', 'url(' + background_image_urls[index] + ')');
                index++;
                image_div.animate({opacity: 1}, 600, function () {
                    if (index == background_image_urls.length) index = 0;
                });
            });
        }, 5000);
    }
});
