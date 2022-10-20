<template>
  <v-container fluid style="max-width: 700px">
    <v-row class="text-left">
      <v-col cols="12">
        <h1 class="secondary--text my-3">Research Interests</h1>
        <p class="grey-darken-1--text">
          Source separation, Independent Component Analysis (ICA), declipping,
          audio inpainting, convex optimisation, probabilistic modelling,
          statistical signal processing, trans and queer marginalisation in
          technology.
        </p>
        <h2 class="secondary--text my-3">Publications</h2>
        <v-list>
          <template v-for="(item, index) in publications">
            <v-list-item :key="index" dark>
              <v-list-item-icon>
                <v-icon color="grey darken-1">chevron_right</v-icon>
              </v-list-item-icon>
              <v-list-item-content>
                <span class="grey-darken-1--text">
                  {{ item.author }}
                  <a
                    :href="item.titleUrl"
                    target="_blank"
                    class="grey-darken-1--link italic"
                  >
                    {{ item.title }}
                  </a>
                  .
                  <a
                    :href="item.conferenceUrl"
                    target="_blank"
                    class="grey-darken-1--link italic"
                  >
                    {{ item.conference }}
                  </a>
                  , {{ item.location }}
                </span>
              </v-list-item-content>
            </v-list-item>
          </template>
        </v-list>
        <p class="grey-darken-1--text">
          (I have changed my legal surname twice. This is why these publications
          show my former names.)
        </p>
        <h2 class="secondary--text my-3">Demonstrating</h2>
        <v-list dense>
          <v-list-item dark>
            <v-list-item-icon>
              <v-icon color="grey darken-1">chevron_right</v-icon>
            </v-list-item-icon>
            <v-list-item-content>
              <span class="grey-darken-1--text">
                2019 to present: None due to leave in order to access
                gender-affirming healthcare
              </span>
            </v-list-item-content>
          </v-list-item>
          <template v-for="(item, index) in demonstrations">
            <v-list-item :key="index" dark>
              <v-list-item-icon>
                <v-icon color="grey darken-1">chevron_right</v-icon>
              </v-list-item-icon>
              <v-list-item-content>
                <span class="grey-darken-1--text">
                  {{ item.title }} (<a
                    :href="item.url"
                    target="_blank"
                    class="grey-darken-1--link italic"
                  >
                    {{ item.courseNo }}
                  </a>
                  )
                </span>
              </v-list-item-content>
            </v-list-item>
          </template>
        </v-list>
        <h2 class="secondary--text my-3">Slides</h2>
        <v-row v-for="(item, index) in slides" :key="index">
          <v-col
            cols="12"
            md="6"
            style="
              display: flex;
              justify-content: center;
              flex-direction: column;
            "
          >
            <span class="grey-darken-1--text">
              <strong> {{ item.title }} </strong>
            </span>
            <span v-if="item.url" class="grey-darken-1--text">
              <strong>
                (<a
                  :href="item.url"
                  target="_blank"
                  class="grey-darken-1--link italic"
                >
                  {{ item.urlText }}
                </a>
                )
              </strong>
              <br />
            </span>
            <span class="grey-darken-1--text">
              {{ item.date }}
            </span>
          </v-col>
          <v-col :key="index" cols="12" md="6">
            <a :href="item.file" target="_blank">
              <v-img
                :src="item.thumbnail"
                height="206px"
                width="275px"
                style="border: black 1px solid"
              />
            </a>
          </v-col>
        </v-row>
      </v-col>
    </v-row>
  </v-container>
</template>
<script>
export default {
  name: "ResearchPanel",
  data() {
    return {
      publications: [
        {
          author: "Charles, A.",
          title:
            "Underdetermined Blind Source Separation for Hard Clipped Stereophonic Mixtures",
          titleUrl: "https://ieeexplore.ieee.org/document/8955550",
          conference:
            "2019 International Conference on Data Mining Workshops (ICDMW)",
          conferenceUrl: "http://icdm2019.bigke.org",
          location: "Beijing, China, 8-11 November 2019",
        },
        {
          author: "Turl, A.",
          title:
            "Joint Blind Source Separation and Declipping: A Geometric Approach for Time Disjoint Sources",
          titleUrl: "geo_isspit.pdf",
          conference:
            "The 17th IEEE International Symposium on Signal Processing and Information Technology (ISSPIT 2017)",
          conferenceUrl: "http://www.isspit.org/isspit/2017/",
          location: "Bilbao, Spain, 18-20 December 2017",
        },
      ],
      demonstrations: [
        {
          courseNo: "26956",
          title: 'Autumn 2018: TA for "MSc Software Workshop"',
          url: "https://www.cs.bham.ac.uk/internal/modules/2018/06-26956/",
        },
        {
          courseNo: "23069",
          title: 'Autumn 2017: Demonstrator for "Introduction to AI"',
          url: "https://www.cs.bham.ac.uk/internal/modules/2017/06-23069/",
        },
        {
          courseNo: "28345",
          title:
            'Autumn 2017: Demonstrator for "Data Structures and Algorithms"',
          url: "https://www.cs.bham.ac.uk/internal/modules/2017/06-28345/",
        },
        {
          courseNo: "23069",
          title: 'Autumn 2016: Demonstrator for "Introduction to AI"',
          url: "https://www.cs.bham.ac.uk/internal/modules/2016/06-23069/",
        },
        {
          courseNo: "28345",
          title:
            'Autumn 2016: Demonstrator for "Data Structures and Algorithms"',
          url: "https://www.cs.bham.ac.uk/internal/modules/2016/06-28345/",
        },
        {
          courseNo: "23069",
          title: 'Autumn 2015: Demonstrator for "Introduction to AI"',
          url: "https://www.cs.bham.ac.uk/internal/modules/2015/06-23069/",
        },
      ],
      slides: [
        {
          date: "18 December 2017",
          file: "isspit17.pdf",
          thumbnail: require("../assets/isspit17.png"),
          title:
            "Joint Blind Source Separation and Declipping: A Geometric Approach for Time Disjoint Sources",
          url: "http://www.isspit.org/isspit/2017/",
          urlText: "ISSPIT 2017",
        },
        {
          date: "13 January 2016",
          file: "ica.pdf",
          thumbnail: require("../assets/ica.jpg"),
          title: "Independent Component Analysis (ICA): An Introduction",
        },
        {
          date: "19 November 2015",
          file: "em_algorithm.pdf",
          thumbnail: require("../assets/em.png"),
          title: "The Expectation-Maximisation (EM) algorithm",
        },
      ],
    };
  },
  methods: {
    goToUrl(url) {
      window.open(url, "_blank");
    },
  },
};
</script>
