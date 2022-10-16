<template>
  <v-app>
    <v-navigation-drawer app right color="secondary" dark>
      <template v-slot:prepend>
        <v-list-item two-line>
          <v-list-item-avatar size="224">
            <img src="./assets/avatar.jpg" />
          </v-list-item-avatar>
        </v-list-item>
        <v-list-item three-line color="primary" dark>
          <v-list-item-content>
            <v-list-item-title style="font-size: 1.25rem"
              >Adahlia S. Lux</v-list-item-title
            >
            <v-list-item-subtitle class="primary--text"
              >Final year PhD student, supervised by
              <a href="http://www.cs.bham.ac.uk/~axk/" target="_blank"
                >Prof. Ata Kab&aacute;n
              </a></v-list-item-subtitle
            >
          </v-list-item-content>
        </v-list-item>
      </template>
      <v-divider color="accent"></v-divider>

      <template v-for="(item, index) in externalLinkItems">
        <v-divider
          v-if="item.divider"
          :key="index"
          color="accent"
          inset
        ></v-divider>
        <v-list-item
          v-else
          :key="index"
          color="primary"
          dark
          @click="goToUrl(item.url)"
        >
          <v-list-item-icon>
            <v-icon>{{ item.icon }}</v-icon>
          </v-list-item-icon>

          <v-list-item-content>
            <v-list-item-title>{{ item.title }}</v-list-item-title>
          </v-list-item-content>
        </v-list-item>
      </template>
    </v-navigation-drawer>

    <v-app-bar flat app dense color="secondary">
      <v-tabs center-active dark v-model="selectedTab">
        <v-spacer />
        <v-tab>About Me</v-tab>
        <v-tab>Research</v-tab>
        <v-tab>Art</v-tab>
        <v-tab>Poetry</v-tab>
        <v-tab>Music</v-tab>
        <v-tab>Contact</v-tab>
      </v-tabs>
    </v-app-bar>

    <v-main>
      <v-tabs-items v-model="selectedTab">
        <v-tab-item><AboutMePanel /></v-tab-item>
        <v-tab-item><ResearchPanel /></v-tab-item>
        <v-tab-item><ArtPanel /></v-tab-item>
        <v-tab-item><PoetryPanel /></v-tab-item>
        <v-tab-item><MusicPanel /></v-tab-item>
        <v-tab-item><ContactPanel /></v-tab-item>
      </v-tabs-items>
    </v-main>
  </v-app>
</template>

<script>
import AboutMePanel from "./components/AboutMePanel";
import ArtPanel from "./components/ArtPanel";
import ContactPanel from "./components/ContactPanel";
import MusicPanel from "./components/MusicPanel";
import PoetryPanel from "./components/PoetryPanel";
import ResearchPanel from "./components/ResearchPanel";

export default {
  name: "App",

  components: {
    AboutMePanel,
    ArtPanel,
    ContactPanel,
    MusicPanel,
    PoetryPanel,
    ResearchPanel,
  },

  metaInfo: {
    title: "Home Page",
    titleTemplate: "%s - Adahlia Lux",
  },

  data: () => ({
    externalLinkItems: [
      {
        icon: "work",
        title: "LinkedIn",
        url: "https://www.linkedin.com/in/adahlia-lux-27121075/",
      },
      { divider: true },
      {
        icon: "smart_display",
        title: "My story on Youtube",
        url: "https://www.youtube.com/watch?v=qEXP7SxQuTs",
      },
      { divider: true },
      {
        icon: "brush",
        title: "My art on Instagram",
        url: "https://www.instagram.com/gotta.lux/",
      },
      { divider: true },
      {
        icon: "coffee",
        title: "Buy me a coffee",
        url: "https://ko-fi.com/gottalux",
      },
      { divider: true },
      {
        icon: "collections",
        title: "Buy my prints",
        url: "https://www.etsy.com/uk/shop/GottaLux",
      },
      { divider: true },
      {
        icon: "shopping_cart",
        title: "Buy my merch",
        url: "https://www.redbubble.com/people/adahlia/shop?asc=u",
      },
      { divider: true },
      {
        icon: "volunteer_activism",
        title: "Fundraiser",
        url: "https://gofund.me/d3c36340",
      },
    ],
    selectedTab: 0,
  }),

  methods: {
    goToUrl(url) {
      window.open(url, "_blank");
    },
  },
};
</script>
<style lang="scss">
@import "./scss/main.scss";
@font-face {
  font-family: "Questrial-Regular";

  src: local("Questrial-Regular"),
    url("./scss/Questrial-Regular.ttf.woff") format("woff");
  font-weight: normal;
  font-style: normal;
}

.external-avatar {
  border-radius: 10%;
  border-color: #fcfcfc;
  border-style: solid;
  border: 3px;
}
</style>
